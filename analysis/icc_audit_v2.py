#!/usr/bin/env python3
"""ICC(1,1) and within-protein mean-square ratio from the cross-evaluation
predictions, for rewriting the ICC paragraph against data.

Naming convention (confirmed from pdisjoint_runs_v2):

    {task}__{train}__{cond}__split{n}[__{source}]__on_{eval}.pred.tsv

    train  : replica | rebuilt      (which construction the model trained on)
    cond   : baseline | ppi | shuffled
    source : absent  -> node2vec PPI embedding
             esm     -> ESM-2 embedding
             prott5  -> ProtT5-XL-U50 embedding
    eval   : replica | rebuilt      (which construction it was scored on)

    8 tasks x 2 train x 2 eval x 3 splits x 7 (cond, source) = 672 files.

`baseline` carries no source and serves as the common sequence-only reference
for all three feature families within the same (task, train, eval, split).

Columns are taken as protein / pos / y / y_pred, matching alldata_runs.  They
are asserted, never inferred from the column name: selecting columns by prefix
is how the site index once ended up averaged in with the predicted scores.

Usage
-----
    python icc_audit_v2.py --root ~/HRP/pdisjoint_runs_v2 --dry-run
    python icc_audit_v2.py --root ~/HRP/pdisjoint_runs_v2 -o icc_audit.tsv
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

PATTERN = (r"^(?P<task>[a-z_0-9]+?)__(?P<train>replica|rebuilt)"
           r"__(?P<cond>baseline|ppi|shuffled)__split(?P<split>\d+)"
           r"(?:__(?P<source>esm|prott5))?"
           r"__on_(?P<eval>replica|rebuilt)\.pred\.tsv$")

PROTEIN, SITE, LABEL, SCORE = "protein", "pos", "y", "y_pred"

ORDER = ["phosphorylation_st", "phosphorylation_y", "acetylation_k",
         "methylation_k", "methylation_r", "sumoylation_k",
         "ubiquitination_k", "glycosylation_n"]


def family(cond, source):
    """(cond, source) -> feature family label.  Explicit table, no prefix or
    suffix matching: `ppi` is the condition flag, not the feature source."""
    if cond == "baseline":
        return "baseline"
    src = {None: "interaction", "esm": "esm2", "prott5": "prott5"}[source]
    return src if cond == "ppi" else src + "_perm"


def icc11(groups):
    """One-way random-effects ICC(1,1) with the Shrout-Fleiss n0 correction for
    unequal group sizes.  Returns (ICC, MSW, MSB, n_groups, n_obs)."""
    k = np.array([len(g) for g in groups], dtype=float)
    a, N = len(groups), k.sum()
    means = np.array([g.mean() for g in groups])
    grand = np.concatenate(groups).mean()
    msb = (k * (means - grand) ** 2).sum() / (a - 1)
    msw = sum(((g - g.mean()) ** 2).sum() for g in groups) / (N - a)
    n0 = (N - (k ** 2).sum() / N) / (a - 1)
    return (msb - msw) / (msb + (n0 - 1) * msw), msw, msb, a, int(N)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.expanduser("~/HRP/pdisjoint_runs_v2"))
    ap.add_argument("--min-sites", type=int, default=5,
                    help="test proteins with fewer sites are dropped (paper uses 5)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-o", "--out", default="icc_audit.tsv")
    args = ap.parse_args()

    paths = sorted(glob.glob(os.path.join(args.root, "*.pred.tsv")))
    parsed, unmatched = [], []
    for p in paths:
        m = re.match(PATTERN, os.path.basename(p))
        if not m:
            unmatched.append(os.path.basename(p))
            continue
        d = m.groupdict()
        d["family"] = family(d["cond"], d["source"])
        d["path"] = p
        parsed.append(d)

    if args.dry_run:
        df = pd.DataFrame(parsed)
        print(f"matched {len(parsed)} / {len(paths)} files")
        if len(df):
            print("\ncounts by (family, train, eval):")
            print(df.groupby(["family", "train", "eval"]).size().to_string())
            print("\ntasks:", sorted(df.task.unique()))
            print("splits:", sorted(df.split.unique()))
        for u in unmatched[:15]:
            print("  UNMATCHED:", u)
        if unmatched:
            print(f"  ... {len(unmatched)} unmatched in total")
        return

    if unmatched:
        raise SystemExit(f"{len(unmatched)} files did not parse, e.g. "
                         f"{unmatched[0]!r}; rerun with --dry-run")

    rows = []
    for d in parsed:
        df = pd.read_csv(d["path"], sep="\t")
        missing = {PROTEIN, SITE, LABEL, SCORE} - set(df.columns)
        assert not missing, f"{d['path']}: missing {missing}"
        groups = [v.to_numpy(float) for _, v in df.groupby(PROTEIN)[SCORE]
                  if len(v) >= args.min_sites]
        if len(groups) < 3:
            continue
        icc, msw, msb, a, N = icc11(groups)
        rows.append({k: d[k] for k in
                     ("task", "train", "eval", "split", "family")}
                    | dict(icc=icc, msw=msw, msb=msb, n_proteins=a, n_sites=N))

    long = pd.DataFrame(rows)
    if long.empty:
        raise SystemExit("nothing scored; check --min-sites")

    key = ["task", "train", "eval", "split"]
    base = (long[long.family == "baseline"]
            .set_index(key)[["icc", "msw"]]
            .rename(columns={"icc": "icc_base", "msw": "msw_base"}))
    out = (long[long.family != "baseline"].set_index(key)
           .join(base, how="inner").reset_index())
    out["icc_fold"] = out.icc / out.icc_base
    out["msw_ratio"] = out.msw_base / out.msw
    out = out[key + ["family", "icc_base", "icc", "icc_fold",
                     "msw_base", "msw", "msw_ratio", "n_proteins", "n_sites"]]
    out.to_csv(args.out, sep="\t", index=False, float_format="%.4f")

    f4 = lambda v: f"{v:.4f}"
    # Ratios are formed AFTER averaging over splits, never by averaging the
    # per-split ratios: ICC(1,1) can be near zero or negative, and a mean of
    # per-split folds then diverges.  Per-split ratios stay in the tsv.
    agg = (out.groupby(["family", "train", "eval", "task"])
              [["icc_base", "icc", "msw_base", "msw"]].mean()
              .reset_index())
    agg["icc_fold"] = np.where(agg.icc_base > 0.01,
                               agg.icc / agg.icc_base, np.nan)
    agg["msw_ratio"] = agg.msw_base / agg.msw
    order = {t: i for i, t in enumerate(ORDER)}
    agg["_task_order"] = agg.task.map(order)
    agg = (agg.sort_values(["family", "train", "eval", "_task_order"])
              .drop(columns="_task_order"))

    for (fam, tr, ev), block in agg.groupby(["family", "train", "eval"],
                                            sort=False):
        print(f"\n=== {fam} | train={tr} | eval={ev} "
              f"(mean over {out.split.nunique()} splits) ===")
        print(block[["task", "icc_base", "icc", "icc_fold", "msw_ratio"]]
              .to_string(index=False, float_format=f4))
        print(f"  icc range {block.icc.min():.4f}-{block.icc.max():.4f} | "
              f"baseline {block.icc_base.min():.4f}-{block.icc_base.max():.4f} | "
              f"fold {block.icc_fold.min():.2f}-{block.icc_fold.max():.2f}"
              f"{' (baseline<=0.01 blanked)' if block.icc_fold.isna().any() else ''} | "
              f"msw_ratio max {block.msw_ratio.max():.2f}")

    print(f"\nwrote {args.out}  ({len(out)} rows)")
    print("\nManuscript claims to check against the blocks above:")
    print("  threshold-trained (train=replica) augmented ICC range   0.23 - 0.87")
    print("  unrestricted-trained (train=rebuilt) augmented ICC      0.25 - 0.35")
    print("  methylation sequence-only ICC                           0.11 - 0.13")
    print("  methylation ICC fold-increase                           1.2 - 1.3")
    print("  acetylation K, ICC 0.043 -> 0.559")
    print("  MSW ratio: unrestricted 0.87-1.00 (PPI), 0.86-1.12 (LM);")
    print("             threshold-trained max 7.9 (PPI) and 7.4 (LM)")


if __name__ == "__main__":
    sys.exit(main())
