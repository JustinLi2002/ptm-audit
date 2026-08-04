#!/usr/bin/env python3
"""2 x 2 x 2 cross-evaluation summary from pdisjoint_runs_v2.

Every trained ensemble is evaluated on BOTH reconstructions' test partitions.
The split file is shared, so the test protein set is identical either way;
only which candidate sites appear as negatives differs. Under a
protein-disjoint split neither test set carries lookup leakage.

Headline result: on rebuilt-test, the real embedding is WORSE than a permuted
vector of the same shape in the five tasks where the shortcut is available
(0/5 tasks, 1/15 partitions), while it is better in all three control tasks
(3/3 tasks, 9/9 partitions). Fisher exact on tasks p = 0.018.
A noise channel can only add capacity; a real feature doing more damage than
noise requires that it carries the information needed to execute the learned
rule, which the natural negative distribution inverts.
"""
import glob
import json
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr

BASE = "/home/FCAM/juli/HRP"
PTMS = ['phosphorylation_y', 'phosphorylation_st', 'ubiquitination_k',
        'sumoylation_k', 'acetylation_k', 'methylation_k', 'methylation_r',
        'glycosylation_n']
# pure-positive share of the rebuilt-evaluated test partitions, mean over seeds
PP = {'phosphorylation_y': 3.5, 'phosphorylation_st': 3.8,
      'ubiquitination_k': 12.6, 'sumoylation_k': 16.5, 'acetylation_k': 20.3,
      'methylation_k': 27.8, 'methylation_r': 31.6, 'glycosylation_n': 35.9}


def load(metric='auroc'):
    R = defaultdict(dict)
    for f in glob.glob(f'{BASE}/pdisjoint_runs_v2/*.json'):
        r = json.load(open(f))
        for t in ('replica', 'rebuilt'):
            R[(r['ptm'], r['dataset'], t)][(r['cond'], r['split_seed'])] = \
                r['tests'][t][metric]
    return R


def cell(R, ptm, train, test):
    d = R[(ptm, train, test)]
    b = [d[('baseline', s)] for s in range(3)]
    pi = [d[('ppi', s)] for s in range(3)]
    pm = [d[('shuffled', s)] for s in range(3)]
    return (b,
            [x - y for x, y in zip(pi, b)],
            [x - y for x, y in zip(pm, b)],
            [x - y for x, y in zip(pi, pm)])


def main():
    for metric in ('auroc', 'auprc'):
        R = load(metric)
        print(f"\n{'#'*100}\n{metric.upper()}\n{'#'*100}")
        for train in ('replica', 'rebuilt'):
            for test in ('replica', 'rebuilt'):
                print(f"\n=== train={train}  test={test} ===")
                print(f"{'PTM':20s} {'tPP%':>6s} {'base':>8s} {'dPPI':>9s} "
                      f"{'dPERM':>9s} {'real-perm':>10s} {'d>0':>5s} {'r>p':>5s}")
                dps, dms, tot_rp = [], [], 0
                for p in PTMS:
                    b, dp, dm, rp = cell(R, p, train, test)
                    nrp = sum(1 for x in rp if x > 0)
                    tot_rp += nrp
                    dps.append(np.mean(dp))
                    dms.append(np.mean(dm))
                    print(f"{p:20s} {PP[p]:5.1f}% {np.mean(b):8.4f} "
                          f"{np.mean(dp):+9.4f} {np.mean(dm):+9.4f} "
                          f"{np.mean(rp):+10.4f} "
                          f"{sum(1 for x in dp if x>0):d}/3 {nrp:d}/3")
                ng = [i for i, p in enumerate(PTMS) if p != 'glycosylation_n']
                r8 = spearmanr([PP[p] for p in PTMS], dps).statistic
                r7 = spearmanr([PP[PTMS[i]] for i in ng],
                               [dps[i] for i in ng]).statistic
                print(f"{'MEAN':20s} {'':6s} {'':8s} {np.mean(dps):+9.4f} "
                      f"{np.mean(dms):+9.4f} {'':10s} "
                      f"{sum(1 for x in dps if x>0)}/8 {tot_rp:2d}/24")
                print(f"  rho(dPPI, purepos): n8={r8:+.3f}  n7(excl N-glyc)={r7:+.3f}")


if __name__ == '__main__':
    main()
