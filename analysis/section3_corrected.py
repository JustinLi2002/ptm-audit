#!/usr/bin/env python3
"""Recompute Section 3 diagnostics without training PTM networks.

All outputs go to a new directory. Existing datasets, predictions and results
are read-only inputs. Protein site sets are deduplicated by (accession, position).
Outer-fold controls never see held-out targets. Main probe comparisons use the
same proteins and folds for both feature families. Percentile ranks are computed
within each task/partition/model before averaging over the three partitions.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
from sklearn.neighbors import NearestNeighbors
from threadpoolctl import threadpool_limits

TASKS = ['phosphorylation_st', 'phosphorylation_y', 'acetylation_k',
         'methylation_k', 'methylation_r', 'sumoylation_k',
         'ubiquitination_k', 'glycosylation_n']
TARGETS = dict(zip(TASKS, ['ST', 'Y', 'K', 'KR', 'R', 'K', 'K', 'N']))
GROUP = dict(zip(TASKS, ['phospho', 'phospho', 'acet', 'methyl', 'methyl',
                        'sumo', 'ubi', 'glyc']))
GROUP_TASKS = {'phospho': ['phosphorylation_st', 'phosphorylation_y'],
               'methyl': ['methylation_k', 'methylation_r'],
               'acet': ['acetylation_k'], 'sumo': ['sumoylation_k'],
               'ubi': ['ubiquitination_k'], 'glyc': ['glycosylation_n']}
GROUP_TARGETS = {'phospho': 'STY', 'methyl': 'KR', 'acet': 'K',
                 'sumo': 'K', 'ubi': 'K', 'glyc': 'N'}
SHORT = dict(zip(TASKS, ['Phospho S/T', 'Phospho Y', 'Acetyl K', 'Meth K/R',
                        'Meth R', 'Sumo K', 'Ubiq K', 'N-Glyc N']))
ALPHAS = np.logspace(-2, 4, 13)


def rho(a, b):
    return float(spearmanr(a, b).statistic)


def read_fasta(path):
    seqs = {}; accession = None
    for line in path.open():
        line = line.strip()
        if line.startswith('>'):
            accession = line.split('|')[1]; seqs[accession] = ''
        elif accession is not None:
            seqs[accession] += line
    return seqs


def counts(sites):
    return sites.groupby('protein').size().astype(int)


def load_data(root, out):
    fasta = root / 'deepmvp/DeepMVP/data/swiss_prot_human_20190214.fasta'
    seqs = read_fasta(fasta)
    frames = {src: {} for src in ('replica', 'rebuilt')}
    positive = {}; audit = []
    for task in TASKS:
        for src in frames:
            df = pd.read_csv(root / src / (task + '_all.tsv'), sep='\t',
                             usecols=['protein', 'pos', 'y'])
            assert not df.duplicated(['protein', 'pos']).any(), (src, task, 'duplicate sites')
            frames[src][task] = df
        pos = frames['rebuilt'][task].query('y == 1')[['protein', 'pos']].drop_duplicates()
        pos2 = frames['replica'][task].query('y == 1')[['protein', 'pos']].drop_duplicates()
        assert set(map(tuple, pos.to_numpy())) == set(map(tuple, pos2.to_numpy())), task
        assert set(pos.protein) <= set(seqs), task
        positive[task] = pos
    merged = {}; sizes = {}
    for group, tasks in GROUP_TASKS.items():
        stacked = pd.concat([positive[t] for t in tasks], ignore_index=True)
        union = stacked.drop_duplicates(['protein', 'pos'])
        merged[group] = counts(union)
        old = counts(stacked)
        audit.append(dict(group=group, concatenated_sites=len(stacked),
                          unique_sites=len(union), repeated_sites=len(stacked)-len(union),
                          proteins_with_changed_depth=int((old != merged[group]).sum())))
        sizes[group] = pd.Series({p: sum(seqs[p].count(aa) for aa in GROUP_TARGETS[group])
                                 for p in merged[group].index})
    overlap = set(map(tuple, positive['methylation_r'].to_numpy())) - set(map(tuple, positive['methylation_k'].to_numpy()))
    assert not overlap, 'Methylation R is not fully contained in the merged K/R release'
    pd.DataFrame(audit).to_csv(out / 'depth_count_audit.tsv', sep='\t', index=False)
    return seqs, frames, positive, merged, sizes


def feature(root, family):
    arr = np.load(root / f'notebooks/protein_features_{family}.npy')
    ids = json.loads((root / f'notebooks/protein_ids_{family}.json').read_text())
    assert len(ids) == len(arr)
    # Match the existing training/analysis lookup: last row wins for repeated
    # accession IDs. Missing accessions cannot match a protein and are excluded.
    return arr, {p: i for i, p in enumerate(ids) if isinstance(p, str) and p}


def ridge_probe(root, out, seqs, merged, sizes):
    feat = {f: feature(root, f) for f in ('ppi', 'esm')}
    fold_rows = []; preds = []
    for group in ('phospho', 'methyl', 'acet', 'ubi', 'sumo', 'glyc'):
        cm = sorted(set(merged[group].index) & set(feat['ppi'][1]) & set(feat['esm'][1]))
        y = np.log1p(merged[group].reindex(cm).to_numpy(float))
        c = np.log1p(sizes[group].reindex(cm).to_numpy(float))[:, None]
        length = np.log1p(np.array([len(seqs[p]) for p in cm], dtype=float))[:, None]
        splits = list(KFold(5, shuffle=True, random_state=0).split(cm))
        for family, (arr, idx) in feat.items():
            X = arr[[idx[p] for p in cm]]
            for fold, (tr, te) in enumerate(splits):
                size_model = LinearRegression().fit(c[tr], y[tr])
                length_model = LinearRegression().fit(length[tr], y[tr])
                size_pred = size_model.predict(c[te])
                train_resid = y[tr] - size_model.predict(c[tr])
                test_resid = y[te] - size_pred
                raw = make_pipeline(StandardScaler(), RidgeCV(alphas=ALPHAS))
                residual = make_pipeline(StandardScaler(), RidgeCV(alphas=ALPHAS))
                raw.fit(X[tr], y[tr]); residual.fit(X[tr], train_resid)
                p_raw = raw.predict(X[te]); p_resid = residual.predict(X[te])
                fold_rows.append(dict(group=group, family=family, fold=fold,
                    n_total=len(cm), n_train=len(tr), n_test=len(te),
                    raw_r2=r2_score(y[te], p_raw), size_r2=r2_score(y[te], size_pred),
                    length_r2=r2_score(y[te], length_model.predict(length[te])),
                    residual_r2=r2_score(test_resid, p_resid),
                    combined_r2=r2_score(y[te], size_pred+p_resid),
                    delta_r2=r2_score(y[te], size_pred+p_resid)-r2_score(y[te], size_pred),
                    raw_alpha=raw[-1].alpha_, residual_alpha=residual[-1].alpha_))
                if group == 'acet':
                    for k,j in enumerate(te):
                        preds.append(dict(protein=cm[j], family=family, fold=fold,
                                          observed=y[j], predicted=p_raw[k]))
            print('RIDGE', group, family, len(cm), flush=True)
    df = pd.DataFrame(fold_rows)
    df.to_csv(out/'ridge_folds.tsv', sep='\t', index=False)
    cols=['raw_r2','size_r2','length_r2','residual_r2','combined_r2','delta_r2']
    summary=df.groupby(['group','family'],sort=False)[cols].mean().reset_index()
    ns=df.groupby(['group','family'],sort=False).n_total.first().reset_index()
    summary=summary.merge(ns,on=['group','family'])
    summary.to_csv(out/'ridge_summary.tsv',sep='\t',index=False)
    pd.DataFrame(preds).to_csv(out/'ridge_acetylation_oof.tsv',sep='\t',index=False)


def pred_path(root, task, train, family, split):
    cond = 'baseline' if family == 'baseline' else 'ppi'
    suffix = '' if family in ('baseline','ppi') else '__'+family
    return root / 'pdisjoint_runs_v2' / f'{task}__{train}__{cond}__split{split}{suffix}__on_rebuilt.pred.tsv'


def score_and_neighbours(root, out, seqs, frames, positive, merged):
    features={f:feature(root,f) for f in ('ppi','esm')}
    rows=[]; knn=[]; support=[]; posrate=[]
    for task in TASKS:
        depth=merged[GROUP[task]]
        for split in range(3):
            sp=pd.read_csv(root/f'pdisjoint/split_seed{split}.csv').set_index('accession')['split']
            b=pd.read_csv(pred_path(root,task,'replica','baseline',split),sep='\t')
            base=b.groupby('protein').y_pred.mean()
            assert sp.reindex(base.index).eq('test').all()
            base_depth=depth.reindex(base.index).fillna(0)
            label_rate=b.groupby('protein').y.mean()
            posrate.append(dict(task=task,split=split,n_proteins=len(base),
                                rho=rho(base_depth,label_rate)))
            rates={}
            for src in ('replica','rebuilt'):
                df=frames[src][task]
                train=df[df.protein.map(sp).eq('train')]
                rates[src]=train.groupby('protein').y.mean()
            for family,(arr,idx) in features.items():
                a=pd.read_csv(pred_path(root,task,'replica',family,split),sep='\t')
                aug=a.groupby('protein').y_pred.mean()
                assert base.index.equals(aug.index)
                rows.append(dict(task=task,family=family,split=split,n_proteins=len(base),
                    baseline=rho(base_depth,base),augmented=rho(base_depth,aug),
                    delta=rho(base_depth,aug)-rho(base_depth,base)))
                eligible={src:set(rates[src].index)&set(idx) for src in rates}
                common=sorted(eligible['replica']&eligible['rebuilt'])
                test=sorted(set(base.index)&set(idx))
                assert set(common).isdisjoint(test)
                nn=NearestNeighbors(n_neighbors=min(25,len(common))).fit(arr[[idx[p] for p in common]])
                _,nb=nn.kneighbors(arr[[idx[p] for p in test]])
                support.append(dict(task=task,family=family,split=split,
                    n_train_thr=len(eligible['replica']),n_train_unr=len(eligible['rebuilt']),
                    n_train_shared=len(common),n_test=len(test)))
                item=dict(task=task,family=family,split=split,n_test=len(test))
                for src in ('replica','rebuilt'):
                    neighbour_rate=rates[src].reindex(common).to_numpy()[nb].mean(axis=1)
                    item['baseline_'+src]=rho(base.reindex(test),neighbour_rate)
                    item['augmented_'+src]=rho(aug.reindex(test),neighbour_rate)
                item['difference']=item['augmented_replica']-item['augmented_rebuilt']
                knn.append(item)
        print('SCORES_KNN',task,flush=True)
    for name,records in [('score_depth',rows),('knn',knn),('knn_support',support),('label_rate_depth',posrate)]:
        df=pd.DataFrame(records);df.to_csv(out/(name+'_splits.tsv'),sep='\t',index=False)
        by=['task','family'] if 'family' in df else ['task']
        df.groupby(by,sort=False).mean(numeric_only=True).reset_index().to_csv(out/(name+'_summary.tsv'),sep='\t',index=False)


def quantiles(values,k=5):
    b=pd.qcut(values,k,labels=False,duplicates='drop')
    if b.isna().all(): b=pd.Series(0,index=values.index)
    return b


def depth_strata(root,out,seqs,positive,merged):
    records=[]; summary=[]; graph=[]
    for task in TASKS:
        # These covariates are fixed independently of model scores and training arm.
        depth=merged[GROUP[task]]
        ps=set()
        for split in range(3):
            b=pd.read_csv(pred_path(root,task,'replica','baseline',split),sep='\t',usecols=['protein'])
            ps.update(b.protein)
        ps=sorted(ps)
        meta=pd.DataFrame(index=ps)
        meta['depth']=depth.reindex(ps).fillna(0)
        meta['size']=[sum(seqs[p].count(a) for a in TARGETS[task]) for p in ps]
        meta['dbin']=quantiles(meta.depth); meta['sbin']=quantiles(meta['size'])
        for family in ('ppi','esm','prott5'):
            for train in ('replica','rebuilt'):
                for split in range(3):
                    b=pd.read_csv(pred_path(root,task,train,'baseline',split),sep='\t')
                    a=pd.read_csv(pred_path(root,task,train,family,split),sep='\t')
                    m=b.merge(a,on=['protein','pos','y'],suffixes=('_base','_aug'),validate='one_to_one')
                    assert len(m)==len(b)==len(a)
                    m=m.join(meta,on='protein')
                    m['change']=m.y_pred_aug.rank(pct=True)-m.y_pred_base.rank(pct=True)
                    for db,g in m.groupby('dbin'):
                        graph.append(dict(task=task,family=family,train=train,split=split,
                            depth_bin=int(db),n_proteins=g.protein.nunique(),n_sites=len(g),
                            depth_median=float(meta.loc[meta.dbin==db,'depth'].median()),
                            pos_rate=g.y.mean(),d_pct=g.change.mean()))
                    for (db,sb),g in m.groupby(['dbin','sbin']):
                        records.append(dict(task=task,family=family,train=train,split=split,
                            depth_bin=int(db),size_bin=int(sb),n_proteins=g.protein.nunique(),
                            n_sites=len(g),depth_median=meta.loc[(meta.dbin==db)&(meta.sbin==sb),'depth'].median(),
                            size_median=meta.loc[(meta.dbin==db)&(meta.sbin==sb),'size'].median(),
                            pos_rate=g.y.mean(),d_pct=g.change.mean(),thin=len(g)<200))
                current=[r for r in graph if (r['task'],r['family'],r['train'])==(task,family,train)]
                df=pd.DataFrame(current)
                for split,g in df.groupby('split'):
                    g=g.sort_values('depth_bin')
                    summary.append(dict(task=task,family=family,train=train,split=int(split),
                        n_bins=len(g),shallow=g.iloc[0].d_pct,deep=g.iloc[-1].d_pct,
                        span=g.iloc[0].d_pct-g.iloc[-1].d_pct,
                        shallow_rate=g.iloc[0].pos_rate,deep_rate=g.iloc[-1].pos_rate))
        print('STRATA',task,flush=True)
    pd.DataFrame(records).to_csv(out/'depth_by_size_splits.tsv',sep='\t',index=False)
    pd.DataFrame(graph).to_csv(out/'depth_plot_splits.tsv',sep='\t',index=False)
    df=pd.DataFrame(summary)
    df.to_csv(out/'depth_endpoints_splits.tsv',sep='\t',index=False)
    df.groupby(['task','family','train'],sort=False).mean(numeric_only=True).reset_index().to_csv(out/'depth_endpoints_summary.tsv',sep='\t',index=False)


def plots(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(8.0,7.1))
    colors={'ppi':'#1f77b4','esm':'#c23b3b'}
    labels={'ppi':'Interaction','esm':'ESM-2'}
    ridge=pd.read_csv(out/'ridge_acetylation_oof.tsv',sep='\t')
    rr=pd.read_csv(out/'ridge_summary.tsv',sep='\t')
    for fam,g in ridge.groupby('family'):
        r=rr.query('group == "acet" and family == @fam').raw_r2.iloc[0]
        axes[0,0].scatter(g.observed,g.predicted,s=2,alpha=.13,color=colors[fam],label=f'{labels[fam]}  CV R2={r:.2f}')
    axes[0,0].plot([0,6],[0,6],':',color='0.4');axes[0,0].set(xlim=(0,6),ylim=(0,6),xlabel='Observed log(1 + merged depth)',ylabel='Out-of-fold prediction',title='a  Depth can be predicted from vectors')
    axes[0,0].legend(fontsize=7,markerscale=4,frameon=False)
    sd=pd.read_csv(out/'score_depth_summary.tsv',sep='\t')
    for i,fam in enumerate(('ppi','esm')):
        for _,r in sd.query('family == @fam').iterrows():
            axes[0,1].plot([0+i*.035,1+i*.035],[r.baseline,r.augmented],'-o',color=colors[fam],alpha=.75,ms=3,lw=.8)
    axes[0,1].axhline(0,color='0.4',ls=':',lw=.7)
    axes[0,1].set(xticks=[0,1],xticklabels=['Sequence only','With protein vector'],ylabel='Spearman rho(score, merged depth)',title='b  Scores shift towards shallow proteins')
    knn=pd.read_csv(out/'knn_summary.tsv',sep='\t').query('family == "ppi"').set_index('task')
    for i,t in enumerate(TASKS):
        r=knn.loc[t];axes[1,0].plot([i,i],[r.augmented_replica,r.augmented_rebuilt],color='0.7')
        axes[1,0].scatter(i,r.augmented_replica,color=colors['ppi'],s=16)
        axes[1,0].scatter(i,r.augmented_rebuilt,color='0.3',s=16)
    axes[1,0].axhline(0,color='0.4',ls=':',lw=.7)
    axes[1,0].set(xticks=range(8),xticklabels=[SHORT[t] for t in TASKS],ylabel='rho(score, neighbour positive fraction)',title='c  Fixed neighbours, two label sources')
    axes[1,0].tick_params(axis='x',labelrotation=55,labelsize=7)
    axes[1,0].legend(handles=[Line2D([],[],ls='',marker='o',color=colors['ppi'],label='Threshold-sampled labels'),Line2D([],[],ls='',marker='o',color='0.3',label='Unrestricted labels')],frameon=False,fontsize=6)
    dep=pd.read_csv(out/'depth_plot_splits.tsv',sep='\t').query('family == "esm"')
    cmap=plt.get_cmap('tab10')
    for i,t in enumerate(TASKS):
        for tr,ls in [('replica','-'),('rebuilt','--')]:
            g=dep.query('task == @t and train == @tr').groupby('depth_bin')[['depth_median','d_pct']].mean()
            axes[1,1].plot(g.depth_median,g.d_pct*100,ls=ls,marker='o',ms=2,color=cmap(i),alpha=1 if tr=='replica' else .55)
    axes[1,1].axhline(0,color='0.4',ls=':',lw=.7)
    axes[1,1].set(xscale='log',xlabel='Merged annotation depth (bin median)',ylabel='Change in mean percentile (points)',title='d  ESM-2 ranking shifts by depth')
    axes[1,1].legend(handles=[Line2D([],[],color='0.3',ls='-',label='Threshold-sampled training'),Line2D([],[],color='0.3',ls='--',label='Unrestricted training')],frameon=False,fontsize=6)
    fig.tight_layout(pad=1.5);fig.savefig(out/'figure2_corrected.png',dpi=240);plt.close(fig)
    strata=pd.read_csv(out/'depth_by_size_splits.tsv',sep='\t').query('family == "esm"')
    fig,axs=plt.subplots(2,8,figsize=(13.0,4.7),sharex=True)
    for i,t in enumerate(TASKS):
        for j,tr in enumerate(('replica','rebuilt')):
            g=strata.query('task == @t and train == @tr').groupby(['depth_bin','size_bin']).d_pct.mean().unstack()
            im=axs[j,i].imshow(g*100,aspect='auto',vmin=-35,vmax=35,cmap='RdBu_r')
            axs[j,i].set_title(SHORT[t] if j==0 else '',fontsize=8)
            axs[j,i].set_yticks(range(len(g)));axs[j,i].set_yticklabels(range(1,len(g)+1),fontsize=6)
            axs[j,i].set_xticks(range(len(g.columns)));axs[j,i].set_xticklabels(range(1,len(g.columns)+1),fontsize=6)
    axs[0,0].set_ylabel('Threshold-trained\nDepth bin',fontsize=8)
    axs[1,0].set_ylabel('Unrestricted-trained\nDepth bin',fontsize=8)
    fig.supxlabel('Full candidate-residue count bin (small to large)',fontsize=9)
    fig.subplots_adjust(left=.065,right=.94,bottom=.16,top=.88,wspace=.3,hspace=.22)
    cax=fig.add_axes([.955,.21,.012,.6]);fig.colorbar(im,cax=cax,label='Percentile-point change')
    fig.savefig(out/'figureS4_corrected.png',dpi=220);plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--skip-ridge',action='store_true')
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    seqs,frames,positive,merged,sizes=load_data(args.root,args.out)
    if not args.skip_ridge: ridge_probe(args.root,args.out,seqs,merged,sizes)
    score_and_neighbours(args.root,args.out,seqs,frames,positive,merged)
    depth_strata(args.root,args.out,seqs,positive,merged)
    plots(args.out)
    manifest={'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'input_root':str(args.root),'status':'complete',
              'depth':'union of unique protein-position annotations in each merged type',
              'size':'full target-residue count from reference FASTA',
              'ridge_sample':'intersection of feature coverage; same proteins and folds',
              'ridge_control':'size fit within each outer training fold only',
              'neighbours':'shared training-protein pool and fixed neighbour IDs for both label sources',
              'percentiles':'within each partition separately; then equal mean across three partitions'}
    (args.out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print('SECTION3_COMPLETE',flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=1): main()
