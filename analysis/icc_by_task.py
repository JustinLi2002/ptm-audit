#!/usr/bin/env python3
"""Grouping of predictions within proteins, all eight tasks, both constructions.

ICC(1,1) of predicted scores across sites of the same protein, over test
proteins carrying at least five sites, with the Shrout-Fleiss correction for
unequal group sizes.

The within-protein mean square (MSW) ratio is the more robust of the two
readings: it has a defined null value of 1.0 and is not bounded by the
baseline ICC. Under corrected sampling the ratio stays at 0.87-1.00 in all
eight tasks; under threshold sampling it reaches 7.9.

The manuscript's ICC figures (0.043 -> 0.559, MSW 0.110 -> 0.031, 662
proteins) are reproduced by pdisjoint_runs/acetylation_k__replica__*__split0,
i.e. protein-disjoint replica-trained split 0 -- NOT the published protocol,
where the section currently sits. The section must be moved and the training
construction stated.
"""
import numpy as np
import pandas as pd

BASE = "/home/FCAM/juli/HRP"
PTMS = ['phosphorylation_st', 'phosphorylation_y', 'ubiquitination_k',
        'sumoylation_k', 'acetylation_k', 'methylation_k', 'methylation_r',
        'glycosylation_n']
PP = {'phosphorylation_y': 3.5, 'phosphorylation_st': 3.8,
      'ubiquitination_k': 12.6, 'sumoylation_k': 16.5, 'acetylation_k': 20.3,
      'methylation_k': 27.8, 'methylation_r': 31.6, 'glycosylation_n': 35.9}


def icc1(d, min_sites=5, min_proteins=20):
    k = d.groupby('protein')['y_pred'].size()
    d = d[d.protein.isin(k[k >= min_sites].index)]
    a = d.protein.nunique()
    if a < min_proteins:
        return np.nan, np.nan, a
    grand = d.y_pred.mean()
    ks = d.groupby('protein')['y_pred'].size().values
    ms = d.groupby('protein')['y_pred'].mean().values
    N = ks.sum()
    MSB = (ks * (ms - grand) ** 2).sum() / (a - 1)
    MSW = d.groupby('protein')['y_pred'] \
           .apply(lambda x: ((x - x.mean()) ** 2).sum()).sum() / (N - a)
    n0 = (N - (ks ** 2).sum() / N) / (a - 1)
    return (MSB - MSW) / (MSB + (n0 - 1) * MSW), MSW, a


def main():
    for tr in ('replica', 'rebuilt'):
        print(f"\n=== trained on {tr}, evaluated on rebuilt-test "
              f"(3 splits averaged) ===")
        print(f"{'PTM':20s} {'tPP%':>6s} {'ICC base':>9s} {'ICC +PPI':>9s} "
              f"{'ICC fold':>9s} {'MSW base':>9s} {'MSW +PPI':>9s} "
              f"{'MSW ratio':>10s} {'n prot':>7s}")
        for p in PTMS:
            acc = {}
            for cond in ('baseline', 'ppi'):
                v = [icc1(pd.read_csv(
                    f'{BASE}/pdisjoint_runs_v2/{p}__{tr}__{cond}__'
                    f'split{s}__on_rebuilt.pred.tsv', sep='\t')) for s in range(3)]
                acc[cond] = (np.nanmean([x[0] for x in v]),
                             np.nanmean([x[1] for x in v]),
                             int(np.mean([x[2] for x in v])))
            ib, wb, ab = acc['baseline']
            iq, wq, _ = acc['ppi']
            print(f"{p:20s} {PP[p]:5.1f}% {ib:9.3f} {iq:9.3f} "
                  f"{iq/ib if ib>0 else float('nan'):8.1f}x {wb:9.4f} "
                  f"{wq:9.4f} {wb/wq if wq>0 else float('nan'):9.2f}x {ab:7d}")


if __name__ == '__main__':
    main()
