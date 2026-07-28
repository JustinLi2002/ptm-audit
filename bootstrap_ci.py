#!/usr/bin/env python3
"""
bootstrap_ci.py — paired bootstrap confidence intervals for Tables 2 and 3.

Table 2 needs intervals on four AUROCs per task (identity and sequence model,
each on the full test set and on the subset with label-homogeneous training
proteins removed). Table 3 needs an interval on the recovery ratio

    (permuted - baseline) / (real - baseline)

which is a ratio of differences and therefore has no closed form.

All quantities for a given task are resampled with the SAME site indices in
every iteration, so the intervals are paired and the differences between
conditions are estimated on matched samples.

Usage:
    python bootstrap_ci.py --n-boot 2000
    python bootstrap_ci.py --table 2 --n-boot 5000
"""

import argparse
import collections
import os

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

BASE = "/home/FCAM/juli/HRP"
RETRAIN = f"{BASE}/retrain"
INF = f"{BASE}/inference_out"

PTMS = ["phosphorylation_st", "phosphorylation_y", "acetylation_k",
        "methylation_k", "methylation_r", "sumoylation_k",
        "ubiquitination_k", "glycosylation_n"]
SHORT = {"phosphorylation_st": "Phospho S/T", "phosphorylation_y": "Phospho Y",
         "acetylation_k": "Acetylation K", "methylation_k": "Methylation K/R",
         "methylation_r": "Methylation R", "sumoylation_k": "Sumoylation K",
         "ubiquitination_k": "Ubiquitination K",
         "glycosylation_n": "N-Glycosylation N"}
# Table 3 covers the six non-phosphorylation tasks
T3_PTMS = ["acetylation_k", "methylation_k", "methylation_r",
           "sumoylation_k", "ubiquitination_k", "glycosylation_n"]


def load_train_stats(ptm):
    """Per-protein (neg, pos) counts in the training partition."""
    per = collections.defaultdict(lambda: [0, 0])
    df = pd.read_csv(f"{RETRAIN}/{ptm}_train.tsv", sep="\t",
                     usecols=["protein", "y"])
    for p, y in zip(df["protein"].values, df["y"].values):
        per[p][int(y)] += 1
    return per


def load_test(ptm, conds):
    """Test labels plus one prediction column per requested condition."""
    base = pd.read_csv(f"{INF}/baseline__{ptm}.tsv", sep="\t")
    out = {"protein": base["protein"].values, "y": base["y"].values.astype(int)}
    for c in conds:
        f = f"{INF}/{c}__{ptm}.tsv"
        if not os.path.exists(f):
            return None
        d = pd.read_csv(f, sep="\t")
        if not np.array_equal(d["protein"].values, out["protein"]):
            raise SystemExit(f"{ptm}/{c}: prediction files are not row-aligned")
        out[c] = d["y_pred"].values
    return out


def safe_auc(y, s):
    return roc_auc_score(y, s) if len(np.unique(y)) > 1 else np.nan


def ci(vals, lo=2.5, hi=97.5):
    v = np.asarray([x for x in vals if np.isfinite(x)])
    if v.size == 0:
        return (np.nan, np.nan)
    return float(np.percentile(v, lo)), float(np.percentile(v, hi))


def fmt(point, lohi, dec=4):
    if not np.isfinite(point):
        return "—"
    return f"{point:.{dec}f} [{lohi[0]:.{dec}f}, {lohi[1]:.{dec}f}]"


def table2(n_boot, seed):
    print("\n=== Table 2 — identity baseline and sequence model, 95% CI ===")
    hdr = (f"{'PTM task':20s} {'identity full':>24s} {'identity restricted':>24s} "
           f"{'sequence full':>24s} {'sequence restricted':>24s}")
    print(hdr)
    print("-" * len(hdr))
    for ptm in PTMS:
        d = load_test(ptm, ["baseline"])
        if d is None:
            print(f"{SHORT[ptm]:20s} missing predictions")
            continue
        per = load_train_stats(ptm)
        rate = {p: pos / (neg + pos) for p, (neg, pos) in per.items()}
        homo = {p for p, (neg, pos) in per.items() if neg == 0 or pos == 0}

        y = d["y"]
        seq = d["baseline"]
        ident = np.array([rate.get(p, 0.5) for p in d["protein"]])
        keep = np.array([p in rate and p not in homo for p in d["protein"]])

        pts = [safe_auc(y, ident), safe_auc(y[keep], ident[keep]),
               safe_auc(y, seq), safe_auc(y[keep], seq[keep])]

        rng = np.random.RandomState(seed)
        n = len(y)
        boots = [[], [], [], []]
        idx_keep = np.flatnonzero(keep)
        for _ in range(n_boot):
            i = rng.randint(0, n, n)                     # full test set
            j = idx_keep[rng.randint(0, idx_keep.size, idx_keep.size)]
            boots[0].append(safe_auc(y[i], ident[i]))
            boots[1].append(safe_auc(y[j], ident[j]))
            boots[2].append(safe_auc(y[i], seq[i]))
            boots[3].append(safe_auc(y[j], seq[j]))

        cells = [fmt(pts[k], ci(boots[k])) for k in range(4)]
        print(f"{SHORT[ptm]:20s} {cells[0]:>24s} {cells[1]:>24s} "
              f"{cells[2]:>24s} {cells[3]:>24s}")

        # the decline, estimated on matched resamples
        dec_id = [b0 - b1 for b0, b1 in zip(boots[0], boots[1])]
        dec_sq = [b2 - b3 for b2, b3 in zip(boots[2], boots[3])]
        print(f"{'':20s}   decline: identity {pts[0]-pts[1]:+.4f} "
              f"[{np.percentile(dec_id,2.5):+.4f}, {np.percentile(dec_id,97.5):+.4f}]"
              f"   sequence {pts[2]-pts[3]:+.4f} "
              f"[{np.percentile(dec_sq,2.5):+.4f}, {np.percentile(dec_sq,97.5):+.4f}]")


def table3(n_boot, seed):
    print("\n=== Table 3 — permutation control, recovery ratio with 95% CI ===")
    hdr = (f"{'PTM task':20s} {'real gain':>22s} {'permuted gain':>22s} "
           f"{'recovery %':>22s}")
    print(hdr)
    print("-" * len(hdr))
    for ptm in T3_PTMS:
        d = load_test(ptm, ["baseline", "kinase", "shuffled"])
        if d is None:
            print(f"{SHORT[ptm]:20s} missing predictions "
                  f"(needs baseline, kinase, shuffled)")
            continue
        y, b, k, s = d["y"], d["baseline"], d["kinase"], d["shuffled"]
        g_real = safe_auc(y, k) - safe_auc(y, b)
        g_perm = safe_auc(y, s) - safe_auc(y, b)
        rec = 100 * g_perm / g_real if abs(g_real) > 1e-9 else np.nan

        rng = np.random.RandomState(seed)
        n = len(y)
        br, bp, brec = [], [], []
        for _ in range(n_boot):
            i = rng.randint(0, n, n)
            ab, ak, as_ = safe_auc(y[i], b[i]), safe_auc(y[i], k[i]), safe_auc(y[i], s[i])
            gr, gp = ak - ab, as_ - ab
            br.append(gr)
            bp.append(gp)
            brec.append(100 * gp / gr if abs(gr) > 1e-9 else np.nan)

        print(f"{SHORT[ptm]:20s} {fmt(g_real, ci(br)):>22s} "
              f"{fmt(g_perm, ci(bp)):>22s} "
              f"{fmt(rec, ci(brec), dec=1):>22s}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--table", choices=["2", "3", "both"], default="both")
    a = ap.parse_args()
    print(f"paired bootstrap, {a.n_boot} resamples, seed {a.seed}")
    if a.table in ("2", "both"):
        table2(a.n_boot, a.seed)
    if a.table in ("3", "both"):
        table3(a.n_boot, a.seed)
    print("\nIntervals are percentile bootstrap over test sites. Resampling is "
          "paired within a task,\nso differences between conditions are estimated "
          "on matched samples. Phosphorylation S/T\nis the largest test set "
          "(36,154 sites) and will dominate runtime.")


if __name__ == "__main__":
    main()
