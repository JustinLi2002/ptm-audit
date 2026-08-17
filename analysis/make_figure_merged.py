#!/usr/bin/env python3
"""Merged cross-evaluation figure: both feature families in one 2x2.

Replaces the separate interaction-embedding and language-model panels. The
comparison that matters is between them -- the interaction channel never
crosses zero while the language model channel does -- and that is invisible
when the two are on different pages.

    python make_figure_merged.py --out figures/

Writes figure4.png.
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
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

BASE = os.environ.get("PTM_AUDIT_BASE", "/home/FCAM/juli/HRP")
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
C_PPI, C_ESM, C_PART = '#1f77b4', '#d62728', '#bbbbbb'
C_MARK = '0.30'          # shape channel only; colour means the family
# Suffixed feature families. 'ppi' carries no suffix and holds the
# baselines, so anything unmatched must fall through to it -- which is
# why the list has to be explicit rather than a single '__esm' test.
SUFFIXED = ('prott5', 'esm')

plt.rcParams.update({'font.size': 8, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.dpi': 300,
                     'savefig.bbox': 'tight'})


def load(metric='auroc'):
    R = defaultdict(dict)
    for f in glob.glob(f'{BASE}/pdisjoint_runs_v2/*.json'):
        stem = os.path.basename(f)[:-5]
        feat = next((k for k in SUFFIXED if stem.endswith('__' + k)), 'ppi')
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


def main(out):
    R = load('auroc')
    order = sorted(PTMS, key=lambda p: PP[p])
    cells = [('replica', 'replica'), ('replica', 'rebuilt'),
             ('rebuilt', 'replica'), ('rebuilt', 'rebuilt')]
    names = {'replica': 'threshold-sampled', 'rebuilt': 'naturally sampled'}

    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.4), sharey=True, sharex=True)
    for ax, (tr, te) in zip(axes.ravel(), cells):
        for i, p in enumerate(order):
            for feat, dx, col in (('ppi', -0.17, C_PPI), ('esm', 0.17, C_ESM)):
                dp, dm = deltas(R, p, tr, te, feat)
                loses = np.mean(dp) < np.mean(dm)
                ax.scatter([i + dx] * 3, dp, s=5, color=C_PART, zorder=2)
                ax.scatter(i + dx, np.mean(dp), s=30, color=col, zorder=3,
                           marker=('v' if loses else 'o'),
                           facecolors=('none' if loses else col),
                           edgecolors=col, linewidths=1.3)
        ax.axhline(0, color='k', lw=0.8, zorder=1)
        for i in range(len(order) - 1):
            ax.axvline(i + 0.5, color='#f2f2f2', lw=0.6, zorder=0)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([SHORT[p] for p in order], rotation=60, ha='right',
                           fontsize=6.5)
        ax.set_title(f'train: {names[tr]}    test: {names[te]}',
                     fontsize=7.5, loc='left')
    axes[0, 0].set_ylabel(r'$\Delta$AUROC')
    axes[1, 0].set_ylabel(r'$\Delta$AUROC')

    fig.legend(handles=[
        Patch(facecolor=C_PPI, edgecolor='none', label='interaction embedding'),
        Patch(facecolor=C_ESM, edgecolor='none', label='frozen ESM-2'),
        Line2D([], [], marker='o', ls='', mfc=C_MARK, mec=C_MARK,
               label='real feature beats its permutation'),
        Line2D([], [], marker='v', ls='', mfc='none', mec=C_MARK,
               label='permutation beats the real feature'),
        Line2D([], [], marker='.', ls='', color=C_PART,
               label='individual partition'),
    ], frameon=False, fontsize=6.5, loc='lower center', ncol=3,
        bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(f'{out}/figure4.png')
    plt.close(fig)
    print('wrote', out + '/figure4.png')

    # numbers for the caption
    for tr, te in cells:
        line = f'  train={tr:8s} test={te:8s}'
        for feat in ('ppi', 'esm'):
            d = [np.mean(deltas(R, p, tr, te, feat)[0]) for p in order]
            rp = sum(1 for p in order
                     if np.mean(deltas(R, p, tr, te, feat)[0])
                     > np.mean(deltas(R, p, tr, te, feat)[1]))
            line += f'  {feat}: {np.mean(d):+.4f} ({sum(1 for x in d if x>0)}/8 pos, {rp}/8 r>p)'
        print(line)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='figures')
    a, _ = ap.parse_known_args()
    os.makedirs(a.out, exist_ok=True)
    main(a.out)
