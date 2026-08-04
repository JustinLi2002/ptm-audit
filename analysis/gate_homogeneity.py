#!/usr/bin/env python3
"""Label homogeneity by construction, partition and split seed.

Gate check for the cross-evaluation: reproduces the manuscript's Table 1
replica column from the reconstruction files, and reports the pure-positive
share of each TEST partition, which determines how the cross-evaluation
result must be read.

Two findings recorded here:
  1. Table 1 reports the PURE-POSITIVE share, not the label-homogeneous share
     defined in Methods. They coincide in seven tasks; in phosphorylation Y
     they differ tenfold (35.8% vs 3.6%) because that task has proteins that
     are homogeneous by chance (all-negative) rather than by construction.
  2. The threshold is not a single parameter. It also changes the class ratio
     (replica 16.7-46.4% positive vs rebuilt 16.7%) because the achievable
     negative pool is itself a function of the threshold.
"""
import collections
import pandas as pd

BASE = "/home/FCAM/juli/HRP"
PTMS = ['acetylation_k', 'glycosylation_n', 'methylation_k', 'methylation_r',
        'phosphorylation_st', 'phosphorylation_y', 'sumoylation_k',
        'ubiquitination_k']
# manuscript Table 1, "Replica (>=10)" column
TABLE1_REPLICA = {'phosphorylation_st': 4.0, 'phosphorylation_y': 3.6,
                  'acetylation_k': 20.9, 'methylation_k': 25.5,
                  'methylation_r': 29.6, 'sumoylation_k': 17.2,
                  'ubiquitination_k': 12.5, 'glycosylation_n': 37.5}


def stats(d):
    """n_sites, positive rate, % sites on homogeneous proteins,
    % sites on pure-positive proteins."""
    per = collections.defaultdict(lambda: [0, 0])
    for p, y in zip(d.protein.values, d.y.values):
        per[p][int(y)] += 1
    homo = sum(n + q for n, q in per.values() if n == 0 or q == 0)
    pure = sum(n + q for n, q in per.values() if n == 0)
    N = len(d)
    return N, d.y.mean(), 100 * homo / N, 100 * pure / N


def main():
    for tag in ('replica', 'rebuilt'):
        print(f"\n{'='*88}\n{tag.upper()}\n{'='*88}")
        print(f"{'PTM':20s} {'part':6s} {'seed':>4s} {'n':>9s} {'pos%':>6s} "
              f"{'homo%':>7s} {'purepos%':>8s} {'unmap':>6s}")
        for nm in PTMS:
            full = pd.read_csv(f'{BASE}/{tag}/{nm}_all.tsv', sep='\t',
                               usecols=['protein', 'y'])
            N, pr, ho, pp = stats(full)
            print(f"{nm:20s} {'ALL':6s} {'-':>4s} {N:9d} {100*pr:5.1f}% "
                  f"{ho:6.1f}% {pp:7.1f}% {'-':>6s}")
            for seed in (0, 1, 2):
                sp = pd.read_csv(f'{BASE}/pdisjoint/split_seed{seed}.csv') \
                       .set_index('accession')['split'].to_dict()
                m = full.protein.map(sp)
                unmap = int(m.isna().sum())
                for part in ('train', 'test'):
                    d = full[m == part]
                    if len(d) == 0:
                        print(f"{'':20s} {part:6s} {seed:4d}  EMPTY")
                        continue
                    N, pr, ho, pp = stats(d)
                    print(f"{'':20s} {part:6s} {seed:4d} {N:9d} {100*pr:5.1f}% "
                          f"{ho:6.1f}% {pp:7.1f}% {unmap:6d}")

    print(f"\n{'='*88}\nGATE: replica homo% against manuscript Table 1\n{'='*88}")
    print(f"{'PTM':20s} {'Table1':>8s} {'ALL':>8s} {'train_s0':>9s} {'d(best)':>9s}")
    for nm in PTMS:
        full = pd.read_csv(f'{BASE}/replica/{nm}_all.tsv', sep='\t',
                           usecols=['protein', 'y'])
        _, _, ho_all, _ = stats(full)
        sp = pd.read_csv(f'{BASE}/pdisjoint/split_seed0.csv') \
               .set_index('accession')['split'].to_dict()
        _, _, ho_tr, _ = stats(full[full.protein.map(sp) == 'train'])
        ref = TABLE1_REPLICA[nm]
        best = min([ho_all, ho_tr], key=lambda v: abs(v - ref))
        print(f"{nm:20s} {ref:7.1f}% {ho_all:7.1f}% {ho_tr:8.1f}% {best-ref:+8.1f}")


if __name__ == '__main__':
    main()
