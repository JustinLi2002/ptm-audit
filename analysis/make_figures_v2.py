#!/usr/bin/env python3
"""Figures 4-7 of the revised manuscript, generated from the run outputs.

Unlike make_figures.py, which draws from hard-coded values, this script reads
pdisjoint_runs_v2/ and the feature matrices directly, so every panel is
reproducible from the archived data.

    python make_figures_v2.py --out figures/

Fig 4  annotation-depth mechanism, three panels
Fig 5  cross-evaluation of the interaction channel, 2x2
Fig 6  cross-evaluation of the frozen language model channel, 2x2
Fig 7  variance structure: within-protein scatter and mean-square ratios
"""
import argparse
import glob
import json
import os
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.neighbors import NearestNeighbors
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

BASE = '/home/FCAM/juli/HRP'
PTMS = ['phosphorylation_y', 'phosphorylation_st', 'ubiquitination_k',
        'sumoylation_k', 'acetylation_k', 'methylation_k', 'methylation_r',
        'glycosylation_n']
SHORT = {'phosphorylation_y': 'Phospho Y', 'phosphorylation_st': 'Phospho S/T',
         'ubiquitination_k': 'Ubiq K', 'sumoylation_k': 'Sumo K',
         'acetylation_k': 'Acetyl K', 'methylation_k': 'Meth K/R',
         'methylation_r': 'Meth R', 'glycosylation_n': 'N-Glyc N'}
PP = {'phosphorylation_y': 3.5, 'phosphorylation_st': 3.8,
      'ubiquitination_k': 12.6, 'sumoylation_k': 16.5, 'acetylation_k': 20.3,
      'methylation_k': 27.8, 'methylation_r': 31.6, 'glycosylation_n': 35.9}
MERGED = {'phosphorylation_y': 'phospho', 'phosphorylation_st': 'phospho',
          'methylation_k': 'methyl', 'methylation_r': 'methyl',
          'ubiquitination_k': 'ubi', 'sumoylation_k': 'sumo',
          'acetylation_k': 'acet', 'glycosylation_n': 'glyc'}
K = 25
C_SEQ, C_PPI, C_ESM, C_PERM = '#4c4c4c', '#1f77b4', '#d62728', '#bbbbbb'

plt.rcParams.update({'font.size': 8, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.dpi': 300,
                     'savefig.bbox': 'tight'})


# ------------------------------------------------------------------ loading
def pred_path(ptm, train, cond, split, feat, test):
    suf = '' if cond == 'baseline' or feat == 'ppi' else f'__{feat}'
    return (f'{BASE}/pdisjoint_runs_v2/{ptm}__{train}__{cond}__'
            f'split{split}{suf}__on_{test}.pred.tsv')


def load_metrics(metric='auroc'):
    R = defaultdict(dict)
    for f in glob.glob(f'{BASE}/pdisjoint_runs_v2/*.json'):
        feat = 'esm' if os.path.basename(f)[:-5].endswith('__esm') else 'ppi'
        r = json.load(open(f))
        for t in ('replica', 'rebuilt'):
            R[(r['ptm'], r['dataset'], t, feat)][(r['cond'], r['split_seed'])] = \
                r['tests'][t][metric]
    return R


def deltas(R, ptm, train, test, feat):
    b = [R[(ptm, train, test, 'ppi')][('baseline', s)] for s in range(3)]
    pi = [R[(ptm, train, test, feat)][('ppi', s)] for s in range(3)]
    pm = [R[(ptm, train, test, feat)][('shuffled', s)] for s in range(3)]
    return ([x - y for x, y in zip(pi, b)], [x - y for x, y in zip(pm, b)])


def depths():
    m = defaultdict(lambda: defaultdict(int))
    c = defaultdict(lambda: defaultdict(int))
    for p in PTMS:
        d = pd.read_csv(f'{BASE}/rebuilt/{p}_all.tsv', sep='\t',
                        usecols=['protein', 'y'])
        for prot, n in d[d.y == 1].groupby('protein').size().items():
            m[MERGED[p]][prot] += n
        for prot, n in d.groupby('protein').size().items():
            c[MERGED[p]][prot] += n
    return m, c


def feature(feat):
    f = np.load(f'{BASE}/notebooks/protein_features_{feat}.npy')
    ids = json.load(open(f'{BASE}/notebooks/protein_ids_{feat}.json'))
    return f, {p: i for i, p in enumerate(ids)}


# ------------------------------------------------------------------ Figure 4
def mechanism_chain_figure(out):
    dm, _ = depths()
    fig, ax = plt.subplots(1, 3, figsize=(7.6, 2.9))
    pipe = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 13)))
    kf = KFold(5, shuffle=True, random_state=0)

    # (a) depth prediction
    for feat, col, lab in [('ppi', C_PPI, 'Interaction'), ('esm', C_ESM, 'ESM-2')]:
        f, idx = feature(feat)
        cm = sorted(set(dm['acet']) & set(idx))
        y = np.log1p([dm['acet'][q] for q in cm])
        p = cross_val_predict(pipe, f[[idx[q] for q in cm]], y, cv=kf)
        r2 = 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        ax[0].scatter(y, p, s=1, alpha=.12, color=col, rasterized=True,
                      label=f'{lab}  $R^2$={r2:.2f}')
    lim = [0, max(ax[0].get_xlim()[1], ax[0].get_ylim()[1])]
    ax[0].plot(lim, lim, ls='--', lw=.6, color='k')
    ax[0].set_xlabel('log annotation depth (observed)')
    ax[0].set_ylabel('predicted')
    ax[0].set_title('a  Feature vs annotation depth', loc='left', fontsize=8)
    ax[0].legend(frameon=False, fontsize=6, markerscale=4)

    # (b) score-depth correlation, baseline vs augmented
    for j, (feat, col) in enumerate([('ppi', C_PPI), ('esm', C_ESM)]):
        for p in PTMS:
            d = dm[MERGED[p]]
            vals = {}
            for cond in ('baseline', 'ppi'):
                v = []
                for s in range(3):
                    g = pd.read_csv(pred_path(p, 'replica', cond, s, feat,
                                              'rebuilt'), sep='\t') \
                        .groupby('protein')['y_pred'].mean().reset_index()
                    g['depth'] = g.protein.map(d).fillna(0)
                    v.append(spearmanr(g.depth, g.y_pred).statistic)
                vals[cond] = np.mean(v)
            x = [0 + j * .06, 1 + j * .06]
            ax[1].plot(x, [vals['baseline'], vals['ppi']], color=col,
                       lw=.8, alpha=.7, marker='o', ms=2.5)
    ax[1].axhline(0, color='k', lw=.6, ls='--')
    ax[1].set_xticks([0, 1])
    ax[1].set_xticklabels(['sequence\nonly', 'with protein-\nlevel channel'])
    ax[1].set_ylabel(r'$\rho$(score, annotation depth)')
    ax[1].set_title('b  The channel inverts it', loc='left', fontsize=8)
    ax[1].legend(handles=[Line2D([], [], color=C_PPI, label='Interaction'),
                          Line2D([], [], color=C_ESM, label='ESM-2')],
                 frameon=False, fontsize=6)

    # (c) neighbourhood positive rate, two label sources
    f, idx = feature('ppi')
    order = sorted(PTMS, key=lambda p: PP[p])
    for i, p in enumerate(order):
        out_ = {}
        for src in ('replica', 'rebuilt'):
            v = []
            for s in range(3):
                sp = pd.read_csv(f'{BASE}/pdisjoint/split_seed{s}.csv') \
                       .set_index('accession')['split'].to_dict()
                tr = pd.read_csv(f'{BASE}/{src}/{p}_all.tsv', sep='\t',
                                 usecols=['protein', 'y'])
                tr = tr[tr.protein.map(sp) == 'train']
                pr = tr.groupby('protein')['y'].mean()
                pr = pr[[q in idx for q in pr.index]]
                nn = NearestNeighbors(n_neighbors=min(K, len(pr))) \
                    .fit(f[[idx[q] for q in pr.index]])
                g = pd.read_csv(pred_path(p, 'replica', 'ppi', s, 'ppi',
                                          'rebuilt'), sep='\t') \
                    .groupby('protein')['y_pred'].mean().reset_index()
                g = g[g.protein.isin(idx)]
                _, nb = nn.kneighbors(f[[idx[q] for q in g.protein]])
                v.append(spearmanr(pr.values[nb].mean(axis=1), g.y_pred).statistic)
            out_[src] = np.mean(v)
        ax[2].plot([i, i], [out_['rebuilt'], out_['replica']], color='#cccccc', lw=.8)
        ax[2].scatter(i, out_['replica'], s=14, color=C_PPI, zorder=3)
        ax[2].scatter(i, out_['rebuilt'], s=14, color=C_SEQ, zorder=3)
    ax[2].axhline(0, color='k', lw=.6, ls='--')
    ax[2].set_xticks(range(8))
    ax[2].set_xticklabels([SHORT[p] for p in order], rotation=60, ha='right', fontsize=6)
    ax[2].set_ylabel(r'$\rho$(score, neighbour pos. rate)', fontsize=7)
    ax[2].set_title('c  Neighbourhood labels', loc='left', fontsize=8)
    ax[2].legend(handles=[Line2D([], [], marker='o', ls='', color=C_PPI,
                                 label='threshold-sampled labels'),
                          Line2D([], [], marker='o', ls='', color=C_SEQ,
                                 label='unrestricted labels')],
                 frameon=False, fontsize=5.5, loc='lower left',
                 bbox_to_anchor=(-0.02, -0.02))
    fig.tight_layout(w_pad=1.8)
    fig.savefig(f'{out}/figure3.png')
    plt.close(fig)


# --------------------------------------------------------------- Figures 5/6
def cross_eval_figure(feat, fname, out, title):
    R = load_metrics('auroc')
    fig, axes = plt.subplots(2, 2, figsize=(6.8, 5.0), sharey=True)
    cells = [('replica', 'replica'), ('replica', 'rebuilt'),
             ('rebuilt', 'replica'), ('rebuilt', 'rebuilt')]
    names = {'replica': 'threshold-sampled', 'rebuilt': 'naturally sampled'}
    order = sorted(PTMS, key=lambda p: PP[p])
    for ax, (tr, te) in zip(axes.ravel(), cells):
        for i, p in enumerate(order):
            dp, dm_ = deltas(R, p, tr, te, feat)
            loses = np.mean(dp) < np.mean(dm_)
            ax.scatter([i] * 3, dp, s=6, color=C_PERM, zorder=2)
            ax.scatter(i, np.mean(dp), s=26, zorder=3,
                       color=(C_ESM if loses else C_PPI),
                       marker=('v' if loses else 'o'))
        ax.axhline(0, color='k', lw=.7)
        ax.set_xticks(range(8))
        ax.set_xticklabels([SHORT[p] for p in order], rotation=60,
                           ha='right', fontsize=6)
        ax.set_title(f'train: {names[tr]}\ntest: {names[te]}', fontsize=7, loc='left')
    axes[0, 0].set_ylabel(r'$\Delta$AUROC')
    axes[1, 0].set_ylabel(r'$\Delta$AUROC')
    fig.suptitle(title, fontsize=9, x=.02, ha='left')
    fig.legend(handles=[Line2D([], [], marker='o', ls='', color=C_PPI,
                               label='real feature beats permutation'),
                        Line2D([], [], marker='v', ls='', color=C_ESM,
                               label='permutation beats real feature'),
                        Line2D([], [], marker='o', ls='', color=C_PERM,
                               label='individual partition')],
               frameon=False, fontsize=6, loc='lower center', ncol=3,
               bbox_to_anchor=(.5, -.02))
    fig.tight_layout(rect=[0, .03, 1, .96])
    fig.savefig(f'{out}/{fname}')
    plt.close(fig)


# ------------------------------------------------------------------ Figure 7
def variance_figure(out):
    fig = plt.figure(figsize=(7.2, 3.1))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.1, 1])

    # (a) within-protein scatter, acetylation K
    axa = fig.add_subplot(gs[0])
    frames = {c: pd.read_csv(pred_path('acetylation_k', 'replica', c, 0,
                                       'ppi', 'rebuilt'), sep='\t')
              for c in ('baseline', 'ppi')}
    k = frames['baseline'].groupby('protein').size()
    prots = k[k >= 5].index
    ranked = frames['ppi'][frames['ppi'].protein.isin(prots)] \
        .groupby('protein')['y_pred'].mean().sort_values().index
    show = list(ranked[::max(1, len(ranked) // 18)])[:18]
    for j, (cond, off) in enumerate([('baseline', 0), ('ppi', 1)]):
        d = frames[cond]
        for i, p in enumerate(show):
            v = d[d.protein == p].y_pred.values
            axa.scatter([i + off * 19] * len(v), v, s=2, alpha=.5,
                        color=(C_SEQ if cond == 'baseline' else C_PPI),
                        rasterized=True)
            axa.plot([i + off * 19 - .35, i + off * 19 + .35],
                     [v.mean()] * 2, color='k', lw=.8)
    axa.set_xticks([9, 28])
    axa.set_xticklabels(['sequence only', 'with interaction channel'])
    axa.set_ylabel('predicted score')
    axa.set_title('a  Acetylation K, 18 test proteins', loc='left', fontsize=8)

    # (b) MSW ratio, all tasks
    axb = fig.add_subplot(gs[1])
    def msw(d, min_sites=5):
        k = d.groupby('protein')['y_pred'].size()
        d = d[d.protein.isin(k[k >= min_sites].index)]
        a, N = d.protein.nunique(), len(d)
        return d.groupby('protein')['y_pred'] \
                .apply(lambda x: ((x - x.mean()) ** 2).sum()).sum() / (N - a)
    order = sorted(PTMS, key=lambda p: PP[p])
    style = [('replica', 'ppi', C_PPI, 'o'), ('replica', 'esm', C_ESM, 'o'),
             ('rebuilt', 'ppi', C_PPI, 'x'), ('rebuilt', 'esm', C_ESM, 'x')]
    for tr, feat, col, mk in style:
        ys = []
        for p in order:
            r = [msw(pd.read_csv(pred_path(p, tr, 'baseline', s, feat,
                                           'rebuilt'), sep='\t')) /
                 msw(pd.read_csv(pred_path(p, tr, 'ppi', s, feat,
                                           'rebuilt'), sep='\t'))
                 for s in range(3)]
            ys.append(np.mean(r))
        axb.plot(range(8), ys, marker=mk, ms=4, lw=.8, color=col,
                 ls=('-' if tr == 'replica' else ':'),
                 label=f'{"threshold" if tr=="replica" else "unrestricted"}-trained, '
                       f'{"interaction" if feat=="ppi" else "ESM-2"}')
    axb.axhline(1.0, color='k', lw=.7, ls='--')
    axb.set_yscale('log')
    axb.set_xticks(range(8))
    axb.set_xticklabels([SHORT[p] for p in order], rotation=60, ha='right', fontsize=6)
    axb.set_ylabel('within-protein mean square ratio')
    axb.set_title('b  Variance collapse, all tasks', loc='left', fontsize=8)
    axb.legend(frameon=False, fontsize=5.5, loc='upper left')
    fig.tight_layout()
    fig.savefig(f'{out}/figure5.png')
    plt.close(fig)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='figures')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    print('Figure 3 (mechanism chain) ...', flush=True); mechanism_chain_figure(a.out)
    # The per-family cross-evaluation panels are superseded by
    # analysis/make_figure_merged.py, which puts both feature families
    # into one 2x2 (manuscript Figure 4).
    print('Figure 5 (variance structure) ...', flush=True); variance_figure(a.out)
    print('done ->', a.out)
