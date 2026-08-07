#!/usr/bin/env python3
"""Cross-evaluation summary: 2x2 of training construction against evaluation
construction, for each protein-constant feature family.

Baselines carry no vector at all and are therefore shared between feature
families; they are stored under the 'ppi' key and read from there.

Headline, interaction embedding: on rebuilt-test the real feature is WORSE than
its own permutation.  Headline, language models: the sign of the contribution
depends on the evaluation construction, not on the encoder.

Usage: crosseval_summary.py [--feat {ppi,esm,prott5}]
"""
import argparse
import glob
import json
import os
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr

BASE = "/home/FCAM/juli/HRP"
PTMS = ['phosphorylation_y', 'phosphorylation_st', 'ubiquitination_k',
        'sumoylation_k', 'acetylation_k', 'methylation_k', 'methylation_r',
        'glycosylation_n']
# suffixed feature families; 'ppi' carries no suffix and holds the baselines
SUFFIXED = ('esm', 'prott5')
FEATS = ('ppi',) + SUFFIXED
LABELS = {'ppi': 'interaction embedding',
          'esm': 'frozen ESM-2 650M',
          'prott5': 'frozen ProtT5-XL-U50'}

# pure-positive share of the rebuilt-evaluated test partitions, mean over seeds
PP = {'phosphorylation_y': 3.5, 'phosphorylation_st': 3.8,
      'ubiquitination_k': 12.6, 'sumoylation_k': 16.5, 'acetylation_k': 20.3,
      'methylation_k': 27.8, 'methylation_r': 31.6, 'glycosylation_n': 35.9}


def feat_of(path):
    """Which feature family a run file belongs to.

    Longest suffix wins, so a family whose name ends in another's would not be
    silently misfiled.
    """
    stem = os.path.basename(path)[:-5]
    for k in sorted(SUFFIXED, key=len, reverse=True):
        if stem.endswith('__' + k):
            return k
    return 'ppi'


def load(metric):
    R = defaultdict(dict)
    for f in glob.glob(f'{BASE}/pdisjoint_runs_v2/*.json'):
        ff = feat_of(f)
        r = json.load(open(f))
        for t in ('replica', 'rebuilt'):
            R[(r['ptm'], r['dataset'], t, ff)][(r['cond'], r['split_seed'])] = \
                r['tests'][t][metric]
    return R


def present(R, feat):
    """Feature families with a complete 8 x 2 x 2 x 3 set of runs."""
    for ptm in PTMS:
        for train in ('replica', 'rebuilt'):
            for test in ('replica', 'rebuilt'):
                for cond in ('ppi', 'shuffled'):
                    for s in range(3):
                        if (cond, s) not in R[(ptm, train, test, feat)]:
                            return False
    return True


def cell(R, ptm, train, test, feat):
    b = [R[(ptm, train, test, 'ppi')][('baseline', s)] for s in range(3)]
    pi = [R[(ptm, train, test, feat)][('ppi', s)] for s in range(3)]
    pm = [R[(ptm, train, test, feat)][('shuffled', s)] for s in range(3)]
    return (b, [x - y for x, y in zip(pi, b)],
            [x - y for x, y in zip(pm, b)], [x - y for x, y in zip(pi, pm)])


def per_family(a):
    for metric in ('auroc', 'auprc'):
        R = load(metric)
        lab = LABELS[a.feat]
        print(f"\n{'#' * 100}\n{metric.upper()}   [{lab}]\n{'#' * 100}")
        for train in ('replica', 'rebuilt'):
            for test in ('replica', 'rebuilt'):
                print(f"\n=== train={train}  test={test} ===")
                print(f"{'PTM':20s} {'tPP%':>6s} {'base':>8s} {'dFeat':>9s} "
                      f"{'dPERM':>9s} {'real-perm':>10s} {'d>0':>5s} {'r>p':>5s}")
                dps, dms, tot = [], [], 0
                for p in PTMS:
                    b, dp, dm, rp = cell(R, p, train, test, a.feat)
                    n = sum(1 for x in rp if x > 0)
                    tot += n
                    dps.append(np.mean(dp))
                    dms.append(np.mean(dm))
                    print(f"{p:20s} {PP[p]:5.1f}% {np.mean(b):8.4f} "
                          f"{np.mean(dp):+9.4f} {np.mean(dm):+9.4f} "
                          f"{np.mean(rp):+10.4f} "
                          f"{sum(1 for x in dp if x > 0):d}/3 {n:d}/3")
                ng = [i for i, p in enumerate(PTMS) if p != 'glycosylation_n']
                r8 = spearmanr([PP[p] for p in PTMS], dps).statistic
                r7 = spearmanr([PP[PTMS[i]] for i in ng],
                               [dps[i] for i in ng]).statistic
                print(f"{'MEAN':20s} {'':6s} {'':8s} {np.mean(dps):+9.4f} "
                      f"{np.mean(dms):+9.4f} {'':10s} "
                      f"{sum(1 for x in dps if x > 0)}/8 {tot:2d}/24")
                print(f"  rho(dFeat, purepos): n8={r8:+.3f}  "
                      f"n7(excl N-glyc)={r7:+.3f}")


def across_families(avail):
    """One row per cell, one column per feature family with complete runs."""
    print(f"\n\n{'#' * 78}\nFeature families compared, mean dChannel per cell"
          f"\n{'#' * 78}")
    for metric in ('auroc', 'auprc'):
        R = load(metric)
        print(f"\n{metric.upper()}")
        head = f"{'train':10s} {'test':10s}" + "".join(
            f" {'d' + f:>10s}" for f in avail)
        if len(avail) > 1:
            head += f" {'spread':>10s}"
        print(head)
        for train in ('replica', 'rebuilt'):
            for test in ('replica', 'rebuilt'):
                m = {f: np.mean([np.mean(cell(R, p, train, test, f)[1])
                                 for p in PTMS]) for f in avail}
                line = f"{train:10s} {test:10s}" + "".join(
                    f" {m[f]:+10.4f}" for f in avail)
                if len(avail) > 1:
                    line += f" {max(m.values()) - min(m.values()):10.4f}"
                print(line)

    # per-task agreement between the two language models, if both are present
    lms = [f for f in avail if f in ('esm', 'prott5')]
    if len(lms) == 2:
        R = load('auroc')
        v = {f: [np.mean(cell(R, p, 'replica', 'rebuilt', f)[1])
                 for p in PTMS] for f in lms}
        r = np.corrcoef(v[lms[0]], v[lms[1]])[0, 1]
        rho = spearmanr(v[lms[0]], v[lms[1]]).statistic
        d = [abs(x - y) for x, y in zip(v[lms[0]], v[lms[1]])]
        print(f"\n{lms[0]} vs {lms[1]} on the harmful cell (AUROC, per task):"
              f"\n  Pearson r={r:+.4f}  Spearman rho={rho:+.4f}"
              f"  median |diff|={np.median(d):.4f}  max={max(d):.4f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--feat', default='ppi', choices=list(FEATS))
    ap.add_argument('--all', action='store_true',
                    help='print the per-family tables for every family found')
    a = ap.parse_args()

    R = load('auroc')
    avail = [f for f in FEATS if present(R, f)]
    missing = [f for f in FEATS if f not in avail]
    if missing:
        print(f"[note] incomplete or absent, excluded: {', '.join(missing)}")
    assert avail, 'no feature family has a complete set of runs'

    for f in (avail if a.all else [a.feat]):
        assert f in avail, f'{f} has an incomplete set of runs'
        a.feat = f
        per_family(a)
    across_families(avail)


if __name__ == '__main__':
    main()
