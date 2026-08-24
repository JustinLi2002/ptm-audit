#!/usr/bin/env python3
"""
summarize_alldata.py — results on the authors' final released dataset.

Reads alldata_runs/*.json and reports, per PTM task, the AUROC and AUPRC of the
sequence-only baseline, the interaction-feature model, and the permuted-feature
control, averaged over validation seeds.

The question this answers: the released protein-level split closes the lookup
channel, but the negative-sampling threshold is unchanged, so 35-97% of training
proteins remain label-homogeneous. Does a protein-level feature channel trained
on that data still transfer to unseen proteins?
"""

import argparse
import glob
import json
import os
from collections import defaultdict

import numpy as np

ORDER = ["phosphorylation_st", "phosphorylation_y", "acetylation_k",
         "methylation_k", "methylation_r", "sumoylation_k",
         "ubiquitination_k", "glycosylation_n"]
SHORT = {"phosphorylation_st": "Phospho S/T", "phosphorylation_y": "Phospho Y",
         "acetylation_k": "Acetylation K", "methylation_k": "Methylation K",
         "methylation_r": "Methylation R", "sumoylation_k": "Sumoylation K",
         "ubiquitination_k": "Ubiquitination K",
         "glycosylation_n": "N-Glycosylation N"}
# share of training sites on label-homogeneous proteins, from the audit
HOMO = {"phosphorylation_st": 4.4, "phosphorylation_y": 4.1,
        "acetylation_k": 24.5, "methylation_k": 9.3, "methylation_r": 8.8,
        "sumoylation_k": 20.4, "ubiquitination_k": 13.5,
        "glycosylation_n": 39.3}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=os.path.expanduser("~/HRP/alldata_runs"))
    a = ap.parse_args()

    R = defaultdict(lambda: defaultdict(list))
    for f in glob.glob(os.path.join(a.dir, "*.json")):
        d = json.load(open(f))
        R[d["ptm"]][d["cond"]].append(d)
    n = sum(len(v) for c in R.values() for v in c.values())
    print(f"{n} runs across {len(R)} tasks\n")

    def m(ptm, cond, key):
        # ddof=1: the sample standard deviation, which is what Supplementary
        # Note S6 and the Table S9 caption quote. numpy's default ddof=0 is the
        # population form and, with only two validation seeds, comes out a
        # factor of sqrt(2) smaller -- enough to make the deposited output
        # disagree with the manuscript in the fourth decimal.
        v = [x[key] for x in R[ptm].get(cond, [])]
        return (np.mean(v), np.std(v, ddof=1), len(v)) if v else (np.nan, np.nan, 0)

    print("=== AUROC ===")
    print(f"{'PTM task':20s} {'homo%':>6s} {'baseline':>16s} {'+PPI':>16s} "
          f"{'permuted':>16s} {'Δ(+PPI)':>9s} {'Δ(perm)':>9s}")
    print("-" * 96)
    gains, pgains = [], []
    for p in ORDER:
        if p not in R:
            continue
        b, bs, _ = m(p, "baseline", "auroc")
        q, qs, _ = m(p, "ppi", "auroc")
        s, ss, _ = m(p, "shuffled", "auroc")
        gains.append(q - b)
        pgains.append(s - b)
        print(f"{SHORT[p]:20s} {HOMO.get(p, float('nan')):5.1f}% "
              f"{b:9.4f}±{bs:.4f} {q:9.4f}±{qs:.4f} {s:9.4f}±{ss:.4f} "
              f"{q-b:+9.4f} {s-b:+9.4f}")
    if gains:
        print(f"\n  mean Δ(+PPI) {np.mean(gains):+.4f}   "
              f"positive in {sum(1 for x in gains if x > 0)}/{len(gains)}")
        print(f"  mean Δ(permuted) {np.mean(pgains):+.4f}   "
              f"positive in {sum(1 for x in pgains if x > 0)}/{len(pgains)}")

    print("\n=== AUPRC ===")
    print(f"{'PTM task':20s} {'pos rate':>9s} {'baseline':>10s} {'+PPI':>10s} "
          f"{'permuted':>10s} {'Δ(+PPI)':>9s} {'Δ(perm)':>9s}")
    print("-" * 82)
    ag, apg = [], []
    for p in ORDER:
        if p not in R:
            continue
        b, _, _ = m(p, "baseline", "auprc")
        q, _, _ = m(p, "ppi", "auprc")
        s, _, _ = m(p, "shuffled", "auprc")
        pos = R[p]["baseline"][0]["pos_test"] if R[p].get("baseline") else float("nan")
        ag.append(q - b)
        apg.append(s - b)
        print(f"{SHORT[p]:20s} {pos:9.3f} {b:10.4f} {q:10.4f} {s:10.4f} "
              f"{q-b:+9.4f} {s-b:+9.4f}")
    if ag:
        print(f"\n  mean Δ(+PPI) {np.mean(ag):+.4f}   "
              f"positive in {sum(1 for x in ag if x > 0)}/{len(ag)}")
        print(f"  mean Δ(permuted) {np.mean(apg):+.4f}   "
              f"positive in {sum(1 for x in apg if x > 0)}/{len(apg)}")

    # does the gain track homogeneity?
    xs = [HOMO[p] for p in ORDER if p in R and HOMO.get(p) is not None]
    ys = [m(p, "ppi", "auroc")[0] - m(p, "baseline", "auroc")[0]
          for p in ORDER if p in R and HOMO.get(p) is not None]
    if len(xs) > 2:
        def rank(v):
            s_ = sorted(range(len(v)), key=lambda i: v[i])
            r = [0.0] * len(v)
            for k, i in enumerate(s_):
                r[i] = k + 1
            return r
        rx, ry = rank(xs), rank(ys)
        mx, my = np.mean(rx), np.mean(ry)
        rho = (sum((i - mx) * (j - my) for i, j in zip(rx, ry)) /
               np.sqrt(sum((i - mx) ** 2 for i in rx) *
                       sum((j - my) ** 2 for j in ry)))
        print(f"\n  homogeneity vs Δ(+PPI): Spearman ρ = {rho:.3f} (n={len(xs)})")


if __name__ == "__main__":
    main()
