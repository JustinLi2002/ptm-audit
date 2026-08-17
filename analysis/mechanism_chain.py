#!/usr/bin/env python3
"""The annotation-depth mechanism, measured end to end.

Link 2  embedding -> annotation depth   ridge CV R2 by merged modification type
Link 3  model score -> annotation depth  rank correlation, baseline vs +PPI
kNN     model score -> neighbourhood positive rate in embedding space
Control the same kNN with neighbour labels taken from BOTH constructions

The control is the decisive one. If the channel tracked real modification
propensity, the correlation with rebuilt-labelled neighbourhoods would be
positive. It is negative in the five high-purepos tasks, i.e. the model
follows a constructed label pattern that runs opposite to biology.

Note on link 2: an earlier version reported R2 = -0.252 +- 0.494. That was an
artifact of KFold(shuffle=False) over task-ordered rows, not a property of the
embedding. Fixed with shuffle=True and StandardScaler.
"""
import argparse
import json
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold, cross_val_score
from sklearn.neighbors import NearestNeighbors
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import os

BASE = os.environ.get("PTM_AUDIT_BASE", "/home/FCAM/juli/HRP")
PTMS = ['phosphorylation_y', 'phosphorylation_st', 'ubiquitination_k',
        'sumoylation_k', 'acetylation_k', 'methylation_k', 'methylation_r',
        'glycosylation_n']
PP = {'phosphorylation_y': 3.5, 'phosphorylation_st': 3.8,
      'ubiquitination_k': 12.6, 'sumoylation_k': 16.5, 'acetylation_k': 20.3,
      'methylation_k': 27.8, 'methylation_r': 31.6, 'glycosylation_n': 35.9}
# the donor threshold is evaluated on the MERGED modification type
MERGED = {'phosphorylation_y': 'phospho', 'phosphorylation_st': 'phospho',
          'methylation_k': 'methyl', 'methylation_r': 'methyl',
          'ubiquitination_k': 'ubi', 'sumoylation_k': 'sumo',
          'acetylation_k': 'acet', 'glycosylation_n': 'glyc'}
K = 25

_ap = argparse.ArgumentParser()
_ap.add_argument('--feat', default='ppi', choices=['ppi', 'esm'])
FEAT = _ap.parse_args().feat
SUF = '' if FEAT == 'ppi' else f'__{FEAT}'
feats = np.load(f'{BASE}/notebooks/protein_features_{FEAT}.npy')
ids = json.load(open(f'{BASE}/notebooks/protein_ids_{FEAT}.json'))
idx = {p: i for i, p in enumerate(ids)}


def depths():
    """Annotation counts from the unrestricted reconstruction, so the depth
    itself does not depend on the threshold being studied."""
    g = defaultdict(int)
    m = defaultdict(lambda: defaultdict(int))
    for p in PTMS:
        d = pd.read_csv(f'{BASE}/rebuilt/{p}_all.tsv', sep='\t',
                        usecols=['protein', 'y'])
        for prot, n in d[d.y == 1].groupby('protein').size().items():
            g[prot] += n
            m[MERGED[p]][prot] += n
    return g, m


def link2(dg, dm):
    print("[link 2] ridge(embedding -> log annotation depth), 5-fold CV R2")
    pipe = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 13)))
    kf = KFold(5, shuffle=True, random_state=0)
    common = sorted(set(dg) & set(idx))
    X = feats[[idx[p] for p in common]]
    y = np.log1p([dg[p] for p in common])
    print(f"  depth vs ||emb||  rho = "
          f"{spearmanr([dg[p] for p in common], np.linalg.norm(X, axis=1)).statistic:+.3f}")
    r2 = cross_val_score(pipe, X, y, cv=kf, scoring='r2')
    print(f"  global          R2 = {r2.mean():+.3f} +- {r2.std():.3f}  (n={len(common)})")
    for mt in ('phospho', 'methyl', 'acet', 'ubi', 'sumo', 'glyc'):
        cm = sorted(set(dm[mt]) & set(idx))
        s = cross_val_score(pipe, feats[[idx[p] for p in cm]],
                            np.log1p([dm[mt][p] for p in cm]), cv=kf, scoring='r2')
        print(f"  {mt:8s}        R2 = {s.mean():+.3f} +- {s.std():.3f}  (n={len(cm)})")


def link3(dm):
    print("\n[link 3] score vs merged-type depth, replica-trained on "
          "rebuilt-test, 3 splits averaged")
    print(f"{'PTM':20s} {'baseline':>10s} {'+PPI':>10s} {'delta':>8s}")
    for p in PTMS:
        d = dm[MERGED[p]]
        row = {}
        for cond in ('baseline', 'ppi'):
            v = []
            for s in range(3):
                f = (f'{BASE}/pdisjoint_runs_v2/{p}__replica__{cond}__'
                     f'split{s}' + ('' if cond == 'baseline' else SUF) + '__on_rebuilt.pred.tsv')
                g = pd.read_csv(f, sep='\t').groupby('protein')['y_pred'] \
                      .mean().reset_index()
                g['depth'] = g.protein.map(d).fillna(0)
                v.append(spearmanr(g.depth, g.y_pred).statistic)
            row[cond] = np.mean(v)
        print(f"{p:20s} {row['baseline']:+10.3f} {row['ppi']:+10.3f} "
              f"{row['ppi']-row['baseline']:+8.3f}")


def knn_control():
    print(f"\n[kNN control] score vs k={K} neighbourhood positive rate, "
          "neighbour labels from both constructions")
    print(f"{'PTM':20s} {'tPP%':>6s} {'nbr=replica':>12s} {'nbr=rebuilt':>12s} "
          f"{'diff':>8s}")
    rows = []
    for p in PTMS:
        out = {}
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
                    .fit(feats[[idx[q] for q in pr.index]])
                d = pd.read_csv(f'{BASE}/pdisjoint_runs_v2/{p}__replica__ppi__'
                                f'split{s}{SUF}__on_rebuilt.pred.tsv', sep='\t')
                g = d.groupby('protein')['y_pred'].mean().reset_index()
                g = g[g.protein.isin(idx)]
                _, nb = nn.kneighbors(feats[[idx[q] for q in g.protein]])
                v.append(spearmanr(pr.values[nb].mean(axis=1), g.y_pred).statistic)
            out[src] = np.mean(v)
        rows.append((PP[p], out['replica'] - out['rebuilt']))
        print(f"{p:20s} {PP[p]:5.1f}% {out['replica']:+12.3f} "
              f"{out['rebuilt']:+12.3f} {out['replica']-out['rebuilt']:+8.3f}")
    ng = [r for r, p in zip(rows, PTMS) if p != 'glycosylation_n']
    print(f"\n  rho(diff, purepos) n8 = "
          f"{spearmanr([r[0] for r in rows], [r[1] for r in rows]).statistic:+.3f}"
          f"   n7 = {spearmanr([r[0] for r in ng], [r[1] for r in ng]).statistic:+.3f}")


if __name__ == '__main__':
    dg, dm = depths()
    link2(dg, dm)
    link3(dm)
    knn_control()
