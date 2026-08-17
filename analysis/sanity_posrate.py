#!/usr/bin/env python3
"""Sanity check for the kNN control in mechanism_chain.py.

If per-protein positive rates in the unrestricted reconstruction were nearly
constant, neighbourhood positive rates computed from rebuilt labels would be
measuring noise and the control would be vacuous. They are not: sd 0.13-0.22,
IQR 0.15-0.32, and positive rate correlates with annotation depth at rho
+0.66 to +0.78 in all eight tasks -- the true biological direction, which the
+PPI channel inverts.

Caveat to state in Methods: part of the depth/positive-rate correlation is
mechanical, since positive rate is npos/(npos + sampled negatives) and the
negative count scales with the number of candidate sites. This does not
affect the argument -- protein-level modification propensity is what the
quantity is meant to express -- but it should be stated rather than left for
a reviewer to raise.
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import os

BASE = os.environ.get("PTM_AUDIT_BASE", "/home/FCAM/juli/HRP")
PTMS = ['phosphorylation_y', 'phosphorylation_st', 'ubiquitination_k',
        'sumoylation_k', 'acetylation_k', 'methylation_k', 'methylation_r',
        'glycosylation_n']


def main():
    print(f"{'PTM':20s} {'posrate mean':>12s} {'sd':>7s} {'IQR':>7s} "
          f"{'rho(posrate,depth)':>20s}")
    for p in PTMS:
        d = pd.read_csv(f'{BASE}/rebuilt/{p}_all.tsv', sep='\t',
                        usecols=['protein', 'y'])
        size = d.groupby('protein').size()
        g = d.groupby('protein').agg(pr=('y', 'mean'), npos=('y', 'sum'))
        g = g[size >= 5]
        q = np.percentile(g.pr, [25, 75])
        print(f"{p:20s} {g.pr.mean():12.3f} {g.pr.std():7.3f} {q[1]-q[0]:7.3f} "
              f"{spearmanr(g.pr, g.npos).statistic:+20.3f}")


if __name__ == '__main__':
    main()
