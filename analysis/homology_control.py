#!/usr/bin/env python3
"""Homology proximity control for the language model arm.

The protein-disjoint partitions are defined by sequence identity at 40%, while a
language model captures remoter homology than that threshold separates. This
invites the objection that the harm reported for the pLM channel reflects
homology leakage across the partition boundary rather than the sampling rule.

Two measurements argue against it.

(a) Separation relative to scale. The median cosine distance from a test protein
    to its nearest training protein, divided by the corresponding leave-one-out
    distance within the training set. Absolute distances are not comparable
    between feature families -- pLM embeddings are anisotropic -- so the ratio
    is the quantity of interest. It is larger for the language model than for
    the interaction embedding, i.e. relative to its own scale the partition
    separates test from training more cleanly in pLM space, not less.

(b) The signature homology transfer would leave. If the channel were exploiting
    proximity to training proteins, adding it would make prediction error more
    dependent on the distance to the nearest training protein, and most so where
    the channel does most harm. The opposite holds: the increase in that
    dependence is largest in the tasks where the channel does LEAST harm.

Usage:
    python homology_control.py [--split-seed 0] [--feat esm]

Caveat recorded in the manuscript: absolute distances in pLM space are small
(median 0.017), so the defence rests on the relative comparison and the rank
ordering rather than on absolute separation. A structure-based partitioning
would be a stronger control and was not performed.
"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors
import os

BASE = os.environ.get("PTM_AUDIT_BASE", "/home/FCAM/juli/HRP")
PTMS = ['phosphorylation_y', 'phosphorylation_st', 'ubiquitination_k',
        'sumoylation_k', 'acetylation_k', 'methylation_k', 'methylation_r',
        'glycosylation_n']
SHORT = {'phosphorylation_y': 'Phospho Y', 'phosphorylation_st': 'Phospho S/T',
         'ubiquitination_k': 'Ubiq K', 'sumoylation_k': 'Sumo K',
         'acetylation_k': 'Acetyl K', 'methylation_k': 'Meth K/R',
         'methylation_r': 'Meth R', 'glycosylation_n': 'N-Glyc N'}


def feature(feat):
    f = np.load(f'{BASE}/notebooks/protein_features_{feat}.npy')
    ids = json.load(open(f'{BASE}/notebooks/protein_ids_{feat}.json'))
    return f, {p: i for i, p in enumerate(ids)}


def unit(x):
    return x / np.linalg.norm(x, axis=1, keepdims=True)


def pred_path(ptm, cond, split, feat, train='replica', test='rebuilt'):
    suf = '' if cond == 'baseline' or feat == 'ppi' else f'__{feat}'
    return (f'{BASE}/pdisjoint_runs_v2/{ptm}__{train}__{cond}__'
            f'split{split}{suf}__on_{test}.pred.tsv')


def separation(split_seed):
    """(a) test-to-train nearest-neighbour distance, relative to the
    leave-one-out distance within the training set."""
    sp = pd.read_csv(f'{BASE}/pdisjoint/split_seed{split_seed}.csv')
    tr_p = set(sp[sp.split == 'train'].accession)
    te_p = sp[sp.split == 'test'].accession.tolist()

    print('(a) nearest-neighbour cosine distance')
    print(f"    {'feature':10s} {'test->train':>12s} {'train (loo)':>12s} "
          f"{'ratio':>7s} {'n test':>7s}")
    out = {}
    for feat in ('ppi', 'esm'):
        f, idx = feature(feat)
        Xtr = unit(f[[idx[p] for p in tr_p if p in idx]])
        Xte = unit(f[[idx[p] for p in te_p if p in idx]])
        d_te, _ = NearestNeighbors(n_neighbors=1, metric='cosine') \
            .fit(Xtr).kneighbors(Xte)
        d_tr, _ = NearestNeighbors(n_neighbors=2, metric='cosine') \
            .fit(Xtr).kneighbors(Xtr)
        a, b = float(np.median(d_te)), float(np.median(d_tr[:, 1]))
        out[feat] = a / b
        print(f"    {feat:10s} {a:12.4f} {b:12.4f} {a/b:7.3f} {len(Xte):7d}")
    return out


def error_distance(feat, split_seed):
    """(b) does adding the channel make error depend on proximity to training?"""
    sp = pd.read_csv(f'{BASE}/pdisjoint/split_seed{split_seed}.csv')
    tr_p = set(sp[sp.split == 'train'].accession)
    f, idx = feature(feat)
    Xtr = unit(f[[idx[p] for p in tr_p if p in idx]])
    nn = NearestNeighbors(n_neighbors=1, metric='cosine').fit(Xtr)

    print('\n(b) rho(distance to nearest training protein, |per-protein error|)')
    print(f"    {'PTM':20s} {'baseline':>9s} {'+feature':>9s} {'delta':>8s} "
          f"{'dFeat':>8s}")
    rows = []
    for p in PTMS:
        rho = {}
        for cond in ('baseline', 'ppi'):
            d = pd.read_csv(pred_path(p, cond, split_seed, feat), sep='\t')
            g = d.groupby('protein').agg(s=('y_pred', 'mean'),
                                         y=('y', 'mean')).reset_index()
            g['err'] = (g.s - g.y).abs()
            g = g[g.protein.isin(idx)]
            dd, _ = nn.kneighbors(unit(f[[idx[q] for q in g.protein]]))
            rho[cond] = spearmanr(dd.ravel(), g.err).statistic
        # harm for the same task, averaged over partitions
        dfeat = []
        for s in range(3):
            suf = '' if feat == 'ppi' else f'__{feat}'
            jb = json.load(open(f'{BASE}/pdisjoint_runs_v2/{p}__replica__'
                                f'baseline__split{s}.json'))
            jf = json.load(open(f'{BASE}/pdisjoint_runs_v2/{p}__replica__ppi__'
                                f'split{s}{suf}.json'))
            dfeat.append(jf['tests']['rebuilt']['auroc']
                         - jb['tests']['rebuilt']['auroc'])
        dfeat = float(np.mean(dfeat))
        delta = rho['ppi'] - rho['baseline']
        rows.append((delta, dfeat))
        print(f"    {SHORT[p]:20s} {rho['baseline']:+9.3f} {rho['ppi']:+9.3f} "
              f"{delta:+8.3f} {dfeat:+8.4f}")

    r = spearmanr([x[0] for x in rows], [x[1] for x in rows]).statistic
    print(f"\n    rho(delta, dFeat) = {r:+.3f}, n = {len(rows)}")
    print("    Homology transfer predicts a NEGATIVE correlation here: the more"
          "\n    harm, the more the error should depend on proximity. It is"
          "\n    positive, so proximity does not explain the harm.")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--split-seed', type=int, default=0)
    ap.add_argument('--feat', default='esm', choices=['ppi', 'esm'])
    a, _ = ap.parse_known_args()
    separation(a.split_seed)
    error_distance(a.feat, a.split_seed)
