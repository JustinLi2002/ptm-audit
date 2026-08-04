#!/usr/bin/env python3
"""Cluster (protein-level) bootstrap for Tables 2, S1 and S2.

Replaces the site-level resampling in bootstrap_ci.py, which assumed site
independence -- contradicted by this paper's own ICC of up to 0.87.

Two changes from the earlier version:
  1. Whole proteins are resampled with replacement, taking all their sites.
  2. The test is on the DIFFERENCE of declines (identity minus sequence)
     rather than on whether two marginal intervals overlap. Overlap of
     marginal intervals is a conservative non-standard heuristic; the paired
     difference is the quantity the claim is about.

Counter-intuitively the marginal intervals WIDEN under clustering while the
interval on the difference NARROWS. Restriction is itself a protein-level
operation, so under protein resampling the full and restricted AUROCs move
together; site resampling decoupled them artificially. This must be stated in
Methods.
"""
import numpy as np
import pandas as pd
from scipy.stats import rankdata

BASE = "/home/FCAM/juli/HRP"
B = 2000
PTMS = ['phosphorylation_st', 'phosphorylation_y', 'acetylation_k',
        'methylation_k', 'methylation_r', 'sumoylation_k',
        'ubiquitination_k', 'glycosylation_n']
SHORT = {'phosphorylation_st': 'Phospho S/T', 'phosphorylation_y': 'Phospho Y',
         'acetylation_k': 'Acetyl K', 'methylation_k': 'Meth K/R',
         'methylation_r': 'Meth R', 'sumoylation_k': 'Sumo K',
         'ubiquitination_k': 'Ubiq K', 'glycosylation_n': 'N-Glyc N'}


def auc(y, s):
    n1 = y.sum()
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(s)
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def ci(v):
    v = np.asarray(v)
    v = v[~np.isnan(v)]
    return np.percentile(v, 2.5), np.percentile(v, 97.5)


def main():
    rng = np.random.default_rng(0)
    rows = []
    for ptm in PTMS:
        tr = pd.read_csv(f'{BASE}/retrain/{ptm}_train.tsv', sep='\t',
                         usecols=['protein', 'y'])
        te = pd.read_csv(f'{BASE}/retrain/{ptm}_test.tsv', sep='\t',
                         usecols=['protein', 'y'])
        pred = pd.read_csv(f'{BASE}/inference_out/baseline__{ptm}.tsv', sep='\t',
                           usecols=['protein', 'y', 'y_pred'])
        assert len(pred) == len(te) and (pred.y.values == te.y.values).all(), ptm

        rate = tr.groupby('protein')['y'].mean()
        homo = set(rate[(rate == 0) | (rate == 1)].index)
        seen = set(rate.index)

        d = te.copy()
        d['seq'] = pred.y_pred.values
        d['ident'] = d.protein.map(rate).fillna(0.5).values
        d['keep'] = [(p in seen) and (p not in homo) for p in d.protein]

        prots = d.protein.unique()
        loc = {p: i for i, p in enumerate(prots)}
        d = d.iloc[np.argsort(d.protein.map(loc).values, kind='stable')] \
             .reset_index(drop=True)
        starts = np.searchsorted(d.protein.map(loc).values, np.arange(len(prots)))
        blocks = [np.arange(s, e) for s, e in
                  zip(starts, np.append(starts[1:], len(d)))]

        y, ident, seq, keep = (d.y.values, d.ident.values,
                               d.seq.values, d.keep.values)
        obs = dict(idf=auc(y, ident), sqf=auc(y, seq),
                   idr=auc(y[keep], ident[keep]), sqr=auc(y[keep], seq[keep]))
        obs['id_drop'] = obs['idf'] - obs['idr']
        obs['sq_drop'] = obs['sqf'] - obs['sqr']
        obs['gap'] = obs['id_drop'] - obs['sq_drop']

        bs = {k: [] for k in obs}
        n = len(blocks)
        for _ in range(B):
            i = np.concatenate([blocks[j] for j in rng.integers(0, n, n)])
            yb, ib, sb, kb = y[i], ident[i], seq[i], keep[i]
            v = dict(idf=auc(yb, ib), sqf=auc(yb, sb),
                     idr=auc(yb[kb], ib[kb]), sqr=auc(yb[kb], sb[kb]))
            v['id_drop'] = v['idf'] - v['idr']
            v['sq_drop'] = v['sqf'] - v['sqr']
            v['gap'] = v['id_drop'] - v['sq_drop']
            for k in bs:
                bs[k].append(v[k])
        rows.append((ptm, obs, {k: ci(v) for k, v in bs.items()},
                     len(prots), int(keep.sum()), len(d)))

    print("=" * 104)
    print("TABLE 2 / S1 - cluster (protein) bootstrap, "
          f"{B} resamples, 95% percentile CI")
    print("=" * 104)
    print(f"{'PTM':13s} {'Identity full':>22s} {'Identity restr':>22s} "
          f"{'Sequence full':>22s} {'Sequence restr':>22s}")
    for ptm, o, c, *_ in rows:
        f = lambda k: f"{o[k]:.4f} [{c[k][0]:.4f},{c[k][1]:.4f}]"
        print(f"{SHORT[ptm]:13s} {f('idf'):>22s} {f('idr'):>22s} "
              f"{f('sqf'):>22s} {f('sqr'):>22s}")

    print("\n" + "=" * 104)
    print("TABLE S2 - declines, and the DIFFERENCE of declines")
    print("=" * 104)
    print(f"{'PTM':13s} {'Identity decline':>24s} {'Sequence decline':>24s} "
          f"{'Gap (id-seq)':>24s} {'gap>0':>6s}")
    nsig = 0
    for ptm, o, c, *_ in rows:
        f = lambda k: f"{o[k]:+.4f} [{c[k][0]:+.4f},{c[k][1]:+.4f}]"
        sig = c['gap'][0] > 0
        nsig += sig
        print(f"{SHORT[ptm]:13s} {f('id_drop'):>24s} {f('sq_drop'):>24s} "
              f"{f('gap'):>24s} {'YES' if sig else 'no':>6s}")
    print(f"\ngap CI excludes zero in {nsig}/8 tasks")
    print(f"mean identity decline {np.mean([o['id_drop'] for _,o,*_ in rows]):+.4f}"
          f"   mean sequence decline {np.mean([o['sq_drop'] for _,o,*_ in rows]):+.4f}")
    print(f"\n{'PTM':13s} {'n proteins':>11s} {'n sites':>9s} {'restricted':>17s}")
    for ptm, o, c, npr, nk, nt in rows:
        print(f"{SHORT[ptm]:13s} {npr:11d} {nt:9d} {nk:10d} ({100*nk/nt:.1f}%)")


if __name__ == '__main__':
    main()
