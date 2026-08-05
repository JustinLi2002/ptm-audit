#!/usr/bin/env python3
"""
unseen_protein_eval.py — a zero-cost preview of the protein-disjoint experiment.

Every test set already contains some sites on proteins that never appear in the
training set (0.6%-10.9%, depending on the PTM). For those sites the label-lookup
channel cannot exist: the model has never seen the protein, so a protein-constant
feature vector carries no memorised label frequency.

If the paper's thesis holds, the gain from protein-level features should shrink
sharply on this subset. This script measures that using the predictions already
written by rerun_inference.py — no retraining required.

Caveat: these subsets are small (~120-260 sites each), so per-PTM AUROCs are
noisy. Pooled effects and the direction of change are what to read.

Usage:
    python unseen_protein_eval.py
"""

import os
from collections import defaultdict

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

BASE = "/home/FCAM/juli/HRP"
RETRAIN = f"{BASE}/retrain"
INF = f"{BASE}/inference_out"
ORDER = ['glycosylation_n', 'methylation_r', 'methylation_k', 'acetylation_k',
         'sumoylation_k', 'ubiquitination_k', 'phosphorylation_st',
         'phosphorylation_y']
SHORT = {'glycosylation_n': 'N-Glyc N', 'methylation_r': 'Meth R',
         'methylation_k': 'Meth K', 'acetylation_k': 'Acetyl K',
         'sumoylation_k': 'Sumo K', 'ubiquitination_k': 'Ubiq K',
         'phosphorylation_st': 'Phospho S/T', 'phosphorylation_y': 'Phospho Y'}


def bootstrap_ci(y, p, n_boot=2000, seed=0):
    """Percentile bootstrap CI for AUROC; returns (lo, hi) or (nan, nan)."""
    rng = np.random.RandomState(seed)
    n = len(y)
    if n < 20 or len(set(y)) < 2:
        return float('nan'), float('nan')
    vals = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        if len(set(y[idx])) < 2:
            continue
        vals.append(roc_auc_score(y[idx], p[idx]))
    if not vals:
        return float('nan'), float('nan')
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main():
    print('=== gain from protein-level features, seen vs unseen proteins ===\n')
    hdr = (f"{'PTM':13s} {'n unseen':>9s} {'pos%':>6s} | "
           f"{'base seen':>10s} {'ppi seen':>9s} {'gain':>8s} | "
           f"{'base uns':>9s} {'ppi uns':>8s} {'gain':>8s} {'95% CI of gain':>22s}")
    print(hdr)
    print('-' * len(hdr))

    pooled = defaultdict(list)
    rows = []
    for ptm in ORDER:
        fb = f'{INF}/baseline__{ptm}.tsv'
        fp = f'{INF}/ppi__{ptm}.tsv'
        if not (os.path.exists(fb) and os.path.exists(fp)):
            print(f'{SHORT[ptm]:13s} missing predictions, skipped')
            continue

        train_prots = set(pd.read_csv(f'{RETRAIN}/{ptm}_train.tsv',
                                      sep='\t', usecols=['protein'])['protein'])
        b = pd.read_csv(fb, sep='\t')
        p = pd.read_csv(fp, sep='\t')
        # the two files come from the same test file in the same order
        assert (b['protein'].values == p['protein'].values).all(), \
            f'{ptm}: prediction files are not row-aligned'

        unseen = ~b['protein'].isin(train_prots)
        y = b['y'].values
        yb, yp = b['y_pred'].values, p['y_pred'].values

        def safe(mask, pred):
            yy = y[mask]
            return roc_auc_score(yy, pred[mask]) if len(set(yy)) > 1 else float('nan')

        seen = ~unseen
        bs, ps_ = safe(seen, yb), safe(seen, yp)
        bu, pu = safe(unseen.values, yb), safe(unseen.values, yp)

        # CI for the gain on the unseen subset, by paired bootstrap
        yu = y[unseen.values]
        ybu, ypu = yb[unseen.values], yp[unseen.values]
        rng = np.random.RandomState(0)
        diffs = []
        for _ in range(2000):
            idx = rng.randint(0, len(yu), len(yu))
            if len(set(yu[idx])) < 2:
                continue
            diffs.append(roc_auc_score(yu[idx], ypu[idx]) -
                         roc_auc_score(yu[idx], ybu[idx]))
        ci = (np.percentile(diffs, 2.5), np.percentile(diffs, 97.5)) if diffs \
            else (float('nan'), float('nan'))

        print(f'{SHORT[ptm]:13s} {unseen.sum():9d} {100 * yu.mean():5.1f}% | '
              f'{bs:10.4f} {ps_:9.4f} {ps_ - bs:+8.4f} | '
              f'{bu:9.4f} {pu:8.4f} {pu - bu:+8.4f} '
              f'{f"[{ci[0]:+.3f}, {ci[1]:+.3f}]":>22s}')

        rows.append((ptm, ps_ - bs, pu - bu, int(unseen.sum())))
        pooled['y'].append(yu)
        pooled['b'].append(ybu)
        pooled['p'].append(ypu)

    if rows:
        gs = np.mean([r[1] for r in rows])
        gu = np.mean([r[2] for r in rows])
        print(f'\n  mean gain on seen proteins   : {gs:+.4f}')
        print(f'  mean gain on unseen proteins : {gu:+.4f}'
              f'   ({100 * gu / gs:.0f}% of the seen-protein gain)')

        # pooled across PTMs — ranks are not comparable across tasks, so pool
        # per-PTM AUROCs rather than raw scores; report as a sanity check only
        print('\n  per-PTM detail (gain seen -> gain unseen):')
        for ptm, g_s, g_u, n in rows:
            print(f'    {SHORT[ptm]:13s} {g_s:+.4f} -> {g_u:+.4f}   (n={n})')

    print('\n  NOTE: unseen-protein subsets are small; treat single-PTM values as\n'
          '  indicative. The CD-HIT protein-disjoint retraining is the definitive\n'
          '  test. This is a preview that costs nothing.')


if __name__ == '__main__':
    main()
