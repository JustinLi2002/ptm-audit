from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'render_deps'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

P=Path(__file__).resolve().parents[1]/'results'
TASKS=['phosphorylation_st','phosphorylation_y','acetylation_k','methylation_k','methylation_r','sumoylation_k','ubiquitination_k','glycosylation_n']
SHORT=['Phospho S/T','Phospho Y','Acetyl K','Meth K/R','Meth R','Sumo K','Ubiq K','N-Glyc N']
colors={'ppi':'#23659a','esm':'#b63c39'}
labels={'ppi':'Interaction','esm':'ESM-2'}
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
fig,axes=plt.subplots(2,2,figsize=(8,7.0))
ridge=pd.read_csv(P/'ridge_acetylation_oof.tsv',sep='\t')
rr=pd.read_csv(P/'ridge_summary.tsv',sep='\t')
for fam,g in ridge.groupby('family'):
    r=rr.query('group == "acet" and family == @fam').raw_r2.iloc[0]
    axes[0,0].scatter(g.observed,g.predicted,s=2,alpha=.15,color=colors[fam],label=f'{labels[fam]}  CV $R^2$ = {r:.2f}')
lo=min(0,float(ridge.predicted.min())-.1)
hi=max(float(ridge.observed.max()),float(ridge.predicted.max()))+.1
axes[0,0].plot([lo,hi],[lo,hi],':',color='.4')
axes[0,0].set(xlim=(lo,hi),ylim=(lo,hi),xlabel='Observed log(1 + merged depth)',ylabel='Out-of-fold prediction',title='a  Predicting depth in acetylation K')
axes[0,0].legend(fontsize=7,markerscale=4,frameon=False)
sd=pd.read_csv(P/'score_depth_summary.tsv',sep='\t')
for i,fam in enumerate(('ppi','esm')):
    for _,r in sd.query('family == @fam').iterrows():
        axes[0,1].plot([0+i*.025,1+i*.025],[r.baseline,r.augmented],'-o',color=colors[fam],alpha=.75,ms=3,lw=.8)
axes[0,1].axhline(0,color='.4',ls=':',lw=.7)
axes[0,1].set(xticks=[0,1],xticklabels=['Sequence only','With protein vector'],ylabel='Spearman ρ (score, merged depth)',title='b  Change in score–depth association')
axes[0,1].legend(handles=[Line2D([],[],color=colors[f],label=labels[f]) for f in ('ppi','esm')],frameon=False,fontsize=7)
knn=pd.read_csv(P/'knn_summary.tsv',sep='\t').query('family == "ppi"').set_index('task')
for i,t in enumerate(TASKS):
    r=knn.loc[t]
    axes[1,0].plot([i,i],[r.augmented_replica,r.augmented_rebuilt],color='.7')
    axes[1,0].scatter(i,r.augmented_replica,color=colors['ppi'],s=16)
    axes[1,0].scatter(i,r.augmented_rebuilt,color='.25',s=16)
axes[1,0].axhline(0,color='.4',ls=':',lw=.7)
axes[1,0].set(xticks=range(8),xticklabels=SHORT,ylabel='ρ (score, neighbour positive fraction)',title='c  Fixed neighbours, two label sources')
axes[1,0].tick_params(axis='x',labelrotation=55,labelsize=7)
axes[1,0].legend(handles=[Line2D([],[],ls='',marker='o',color=colors['ppi'],label='Threshold-sampled labels'),Line2D([],[],ls='',marker='o',color='.25',label='Unrestricted labels')],frameon=False,fontsize=6.5,loc='lower left')
dep=pd.read_csv(P/'depth_plot_splits.tsv',sep='\t').query('family == "esm" and task == "acetylation_k"')
for tr,ls,col,lab in [('replica','-',colors['esm'],'Threshold-sampled training'),('rebuilt','--',colors['ppi'],'Unrestricted training')]:
    g=dep.query('train == @tr').groupby('depth_bin')[['depth_median','d_pct']].mean()
    axes[1,1].plot(g.depth_median,g.d_pct*100,ls=ls,marker='o',ms=4,color=col,label=lab)
axes[1,1].axhline(0,color='.4',ls=':',lw=.7)
axes[1,1].set(xscale='log',xlabel='Merged annotation depth (bin median)',ylabel='Change in mean percentile (points)',title='d  ESM-2 ranking shifts in acetylation K')
axes[1,1].legend(frameon=False,fontsize=6.5,loc='lower left')
fig.tight_layout(pad=1.2)
fig.savefig(P/'figure2_corrected.png',dpi=300)
plt.close(fig)

data=pd.read_csv(P/'depth_by_size_splits.tsv',sep='\t').query('family == "esm"')
agg=data.groupby(['task','train','depth_bin','size_bin']).d_pct.mean()*100
lim=np.ceil(agg.abs().max()/5)*5
fig,axs=plt.subplots(4,4,figsize=(7.8,7.6))
for i,(t,name) in enumerate(zip(TASKS,SHORT)):
    for j,(tr,title) in enumerate([('replica','Threshold-trained'),('rebuilt','Unrestricted-trained')]):
        ax=axs[i//2,(i%2)*2+j]
        g=agg.loc[t,tr].unstack()
        im=ax.imshow(g,aspect='auto',vmin=-lim,vmax=lim,cmap='RdBu_r')
        ax.set_title(name+'\n'+title,fontsize=8,pad=4)
        ax.set_yticks(range(len(g)),labels=range(1,len(g)+1),fontsize=7)
        ax.set_xticks(range(len(g.columns)),labels=range(1,len(g.columns)+1),fontsize=7)
        if j==0: ax.set_ylabel('Depth bin',fontsize=8)
        if i//2==3: ax.set_xlabel('Candidate-count bin',fontsize=8)
fig.subplots_adjust(left=.075,right=.985,bottom=.13,top=.945,wspace=.36,hspace=.52)
cax=fig.add_axes([.28,.047,.46,.019])
fig.colorbar(im,cax=cax,orientation='horizontal',label='Change in mean percentile (points)')
fig.savefig(P/'figureS4_corrected.png',dpi=300,bbox_inches='tight',pad_inches=.12)
plt.close(fig)
print('Figures rebuilt. Full heatmap range:',float(agg.min()),float(agg.max()))
