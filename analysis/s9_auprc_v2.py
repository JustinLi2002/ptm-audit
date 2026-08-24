#!/usr/bin/env python3
"""Recompute Table S9 (authors' current release) in both metrics.

v2 fixes a bug in v1: score columns were selected by `startswith("s")`, which
also matched the `site` column, so the site index was averaged in with the
predicted scores.  Score columns are now tracked in an explicit list and never
inferred from the name.

File layout (confirmed):
    ~/HRP/alldata_runs/{task}__{cond}__val{n}__f90.pred.tsv
    columns: protein  aa  pos  y  y_pred
    cond in {baseline, ppi, shuffled}; n in {0, 1}

Each pred.tsv is already the two-initialisation ensemble (its AUROC equals the
`auroc` field of the sibling .json; the `seed_aurocs` field holds the two
individual initialisations).  This script therefore aggregates only over the
two validation seeds, and reports both ways of doing so:

    mean_of_metrics  -- score each val seed, then average the metric
    ensemble         -- average the two predicted scores, then score once

They are not equivalent.  Compare mean_of_metrics against Table S9 to confirm
which convention S9 used before trusting the AUPRC numbers.

Usage
-----
    python s9_auprc_v2.py --root ~/HRP/alldata_runs
    python s9_auprc_v2.py --root ~/HRP/alldata_runs -o s9_both_metrics.tsv
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

PATTERN = (r"^(?P<task>[a-z_]+)__(?P<cond>baseline|ppi|shuffled)"
           r"__val(?P<seed>\d+)__f90\.pred\.tsv$")

PROTEIN, SITE, LABEL, SCORE = "protein", "pos", "y", "y_pred"

ORDER = ["phosphorylation_st", "phosphorylation_y", "acetylation_k",
         "methylation_k", "methylation_r", "sumoylation_k",
         "ubiquitination_k", "glycosylation_n"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.expanduser("~/HRP/alldata_runs"))
    ap.add_argument("-o", "--out", default="s9_both_metrics.tsv")
    args = ap.parse_args()

    cells = {}
    for p in sorted(glob.glob(os.path.join(args.root, "*.pred.tsv"))):
        m = re.match(PATTERN, os.path.basename(p))
        if m:
            d = m.groupdict()
            cells.setdefault((d["task"], d["cond"]), []).append((d["seed"], p))
    if not cells:
        raise SystemExit(f"nothing matched under {args.root!r}")

    rows = []
    for (task, cond), items in sorted(cells.items()):
        items.sort()
        per_seed, frames, score_cols = [], [], []
        for seed, path in items:
            df = pd.read_csv(path, sep="\t")
            missing = {PROTEIN, SITE, LABEL, SCORE} - set(df.columns)
            assert not missing, f"{path}: missing columns {missing}"
            y, s = df[LABEL].to_numpy(int), df[SCORE].to_numpy(float)
            per_seed.append((roc_auc_score(y, s), average_precision_score(y, s)))
            col = f"score_val{seed}"
            score_cols.append(col)
            frames.append(df[[PROTEIN, SITE, LABEL, SCORE]]
                          .rename(columns={SCORE: col}))

        merged = frames[0]
        for f in frames[1:]:
            merged = merged.merge(f, on=[PROTEIN, SITE, LABEL], how="inner",
                                  validate="one_to_one")
        assert len(merged) == len(frames[0]), (
            f"{task}/{cond}: merge lost rows ({len(frames[0])} -> {len(merged)})")

        y = merged[LABEL].to_numpy(int)
        ens = merged[score_cols].to_numpy(float).mean(axis=1)   # explicit list
        a = np.array(per_seed)
        rows.append(dict(
            task=task, cond=cond, n_val_seeds=len(items), n_sites=len(y),
            pos_rate=y.mean(),
            auroc_meanmetric=a[:, 0].mean(), auprc_meanmetric=a[:, 1].mean(),
            auroc_ensemble=roc_auc_score(y, ens),
            auprc_ensemble=average_precision_score(y, ens),
            auroc_sd=a[:, 0].std(ddof=1) if len(a) > 1 else np.nan,
        ))

    long = pd.DataFrame(rows)
    out = {}
    for agg in ("meanmetric", "ensemble"):
        w = long.pivot(index="task", columns="cond",
                       values=[f"auroc_{agg}", f"auprc_{agg}"])
        w.columns = [f"{m.split('_')[0]}_{c}" for m, c in w.columns]
        for met in ("auroc", "auprc"):
            w[f"d_ppi_{met}"] = w[f"{met}_ppi"] - w[f"{met}_baseline"]
            w[f"d_perm_{met}"] = w[f"{met}_shuffled"] - w[f"{met}_baseline"]
            w[f"real_minus_perm_{met}"] = w[f"{met}_ppi"] - w[f"{met}_shuffled"]
        out[agg] = w.reindex([t for t in ORDER if t in w.index])

    fmt = lambda v: f"{v:+.4f}"
    for agg, w in out.items():
        print(f"\n{'='*72}\n{agg}\n{'='*72}")
        print("absolute:")
        print(w[[c for c in w.columns if not c.startswith(("d_", "real_"))]]
              .to_string(float_format=lambda v: f"{v:.4f}"))
        show = [c for c in w.columns if c.startswith(("d_", "real_minus_"))]
        print("\ndeltas:")
        print(w[show].to_string(float_format=fmt))
        print("\nmean:")
        print(w[show].mean().to_string(float_format=fmt))
        print("\nsign counts (of %d tasks):" % len(w))
        for c in show:
            neg = int((w[c] < 0).sum())
            print(f"  {c:26s} negative {neg}, positive {len(w)-neg}")
        print("\nAUPRC vs AUROC sign agreement:")
        for base in ("d_ppi", "d_perm", "real_minus_perm"):
            same = np.sign(w[f"{base}_auroc"]) == np.sign(w[f"{base}_auprc"])
            n = int(same.sum())
            print(f"  {base:16s} {n}/{len(w)}"
                  + ("" if n == len(w) else
                     "   DISAGREES: " + ", ".join(w.index[~same])))

    print(f"\n{'='*72}")
    print("across-val-seed SD of AUROC, by task and condition:")
    print(long.pivot(index="task", columns="cond", values="auroc_sd")
              .reindex([t for t in ORDER if t in set(long.task)])
              .to_string(float_format=lambda v: f"{v:.4f}"))

    out["meanmetric"].to_csv(args.out, sep="\t", float_format="%.4f")
    print(f"\nwrote {args.out}  (mean-of-metrics convention)")


if __name__ == "__main__":
    sys.exit(main())
