#!/usr/bin/env python3
"""Within-task dose-response: does the harm fall on shallowly annotated proteins?

The cross-task rank correlations rest on seven or eight points, and the eight
tasks are not eight independent units. This replaces that argument with a
within-task one: bin the test proteins of a single task by annotation depth and
ask whether the harm concentrates in the shallow bins. The unit becomes
thousands of sites per task rather than eight tasks.

The confound that has to be controlled is protein size. Annotation depth is
strongly correlated with the number of candidate sites -- log candidate count
alone explains 0.132-0.322 of the embedding's variance in the manuscript -- so
binning by depth also bins by size, and a gradient along depth could be a
gradient along size. The default here is therefore a two-way stratification:
depth quantiles crossed with candidate-count quantiles. A real depth effect
appears as a gradient down the columns that persists within every row.

Annotation depth is taken from the unrestricted reconstruction, so it is the
number of sites a protein actually carries, not the number the threshold-sampled
partition happens to expose.

The primary quantity is NOT a within-bin AUPRC. The shortcut acts between
proteins, not within them: a model that has learned "shallow protein implies
positive" raises every site of a shallow protein by roughly the same amount, so
inside a bin of uniformly shallow proteins the ranking barely moves and the
within-bin AUPRC is near zero regardless of how much damage the shortcut does.
Binning by depth and scoring inside the bin removes the very effect it is meant
to measure.

What the shortcut does predict is a shift in where a bin's sites sit in the
ranking of the whole test set. So the primary measure is the change in mean
score percentile of a bin's sites, computed against the full task, alongside
that bin's true positive rate. The prediction is a mismatch: shallow bins gain
percentile while carrying a lower positive rate than the deep bins that lose it.
Within-bin AUPRC is reported as well, and is expected to be weak for the reason
above.

Usage
-----
    python depth_stratified.py --root ~/HRP/pdisjoint_runs_v2 \\
        --data ~/HRP/rebuilt --dry-run
    python depth_stratified.py --root ~/HRP/pdisjoint_runs_v2 \\
        --data ~/HRP/rebuilt --family esm2 -o depth_strata.tsv
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

PATTERN = (r"^(?P<task>[a-z_0-9]+?)__(?P<train>replica|rebuilt)"
           r"__(?P<cond>baseline|ppi|shuffled)__split(?P<split>\d+)"
           r"(?:__(?P<source>esm|prott5))?"
           r"__on_(?P<eval>replica|rebuilt)\.pred\.tsv$")

PROTEIN, SITE, LABEL, SCORE = "protein", "pos", "y", "y_pred"

ORDER = ["phosphorylation_st", "phosphorylation_y", "acetylation_k",
         "methylation_k", "methylation_r", "sumoylation_k",
         "ubiquitination_k", "glycosylation_n"]


def family(cond, source):
    if cond == "baseline":
        return "baseline"
    src = {None: "interaction", "esm": "esm2", "prott5": "prott5"}[source]
    return src if cond == "ppi" else src + "_perm"


def depth_and_size(data_dir, task):
    """Annotation depth and candidate-site count per protein, from the
    unrestricted reconstruction: depth is the number of positive sites the
    protein carries there, size the number of candidate residues."""
    path = os.path.join(data_dir, f"{task}_all.tsv")
    df = pd.read_csv(path, sep="\t")
    cols = {c.lower(): c for c in df.columns}
    p = cols.get("protein") or cols.get("prot")
    y = cols.get("y") or cols.get("label")
    assert p and y, f"{path}: columns are {list(df.columns)}"
    g = df.groupby(df[p])
    return pd.DataFrame({"depth": g[y].sum(), "size": g.size()})


def qbin(v, k, label):
    """Quantile bins that tolerate heavy ties, which depth has: many proteins
    carry exactly one or two sites, so equal-frequency edges collapse."""
    try:
        b = pd.qcut(v, k, labels=False, duplicates="drop")
    except ValueError:
        b = pd.Series(np.zeros(len(v), dtype=int), index=v.index)
    n = int(b.max()) + 1 if len(b) else 0
    if n < k:
        print(f"    note: {label} collapsed to {n} bins (ties)", file=sys.stderr)
    return b


def score(y, s):
    if y.min() == y.max():
        return np.nan, np.nan
    return roc_auc_score(y, s), average_precision_score(y, s)


def percentiles(v):
    """Rank each score as a percentile of the whole task, so that bins can be
    compared on where their sites sit in the ranking the model actually
    produces."""
    return pd.Series(v).rank(pct=True).to_numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.expanduser("~/HRP/pdisjoint_runs_v2"))
    ap.add_argument("--data", default=os.path.expanduser("~/HRP/rebuilt"),
                    help="unrestricted reconstruction, for depth and size")
    ap.add_argument("--family", default="esm2",
                    choices=["interaction", "esm2", "prott5"])
    ap.add_argument("--train", default="replica", choices=["replica", "rebuilt"],
                    help="replica is the threshold-trained arm; rebuilt is the "
                         "negative control, where no gradient should appear")
    ap.add_argument("--eval", default="rebuilt", choices=["replica", "rebuilt"])
    ap.add_argument("--depth-bins", type=int, default=5)
    ap.add_argument("--size-bins", type=int, default=5,
                    help="set to 1 to disable the size stratification")
    ap.add_argument("--min-sites", type=int, default=200,
                    help="bins with fewer pooled sites are reported but flagged")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-o", "--out", default="depth_strata.tsv")
    args = ap.parse_args()

    files = {}
    for p in sorted(glob.glob(os.path.join(args.root, "*.pred.tsv"))):
        m = re.match(PATTERN, os.path.basename(p))
        if not m:
            raise SystemExit(f"did not parse: {os.path.basename(p)}")
        d = m.groupdict()
        if d["train"] == args.train and d["eval"] == args.eval:
            key = (d["task"], family(d["cond"], d["source"]))
            files.setdefault(key, {})[d["split"]] = p

    tasks = [t for t in ORDER if (t, args.family) in files]
    if args.dry_run:
        print(f"family={args.family} train={args.train} eval={args.eval}")
        for t in tasks:
            print(f"  {t:20s} splits={sorted(files[(t, args.family)])}")
        print(f"\ndepth/size source: {args.data}")
        for t in tasks[:2]:
            ds = depth_and_size(args.data, t)
            print(f"  {t}: {len(ds)} proteins, depth "
                  f"{ds.depth.min()}-{ds.depth.max()} (median {ds.depth.median():.0f}), "
                  f"size {ds['size'].min()}-{ds['size'].max()}")
        return

    rows = []
    for task in tasks:
        ds = depth_and_size(args.data, task)
        frames = []
        for s in sorted(files[(task, args.family)]):
            base = pd.read_csv(files[(task, "baseline")][s], sep="\t")
            aug = pd.read_csv(files[(task, args.family)][s], sep="\t")
            m = base.merge(aug, on=[PROTEIN, SITE, LABEL],
                           suffixes=("_base", "_aug"), validate="one_to_one")
            assert len(m) == len(base), f"{task}/{s}: merge lost rows"
            frames.append(m)
        m = pd.concat(frames, ignore_index=True)
        m = m.join(ds, on=PROTEIN)
        missing = m.depth.isna().sum()
        if missing:
            print(f"  {task}: {missing} sites on proteins absent from {args.data}, "
                  f"dropped", file=sys.stderr)
            m = m.dropna(subset=["depth"])

        prot = m.drop_duplicates(PROTEIN).set_index(PROTEIN)[["depth", "size"]]
        prot["dbin"] = qbin(prot.depth, args.depth_bins, f"{task} depth")
        prot["sbin"] = (qbin(prot["size"], args.size_bins, f"{task} size")
                        if args.size_bins > 1 else 0)
        m = m.join(prot[["dbin", "sbin"]], on=PROTEIN)
        # percentiles over the whole task, before any binning
        m["pct_base"] = percentiles(m[f"{SCORE}_base"].to_numpy(float))
        m["pct_aug"] = percentiles(m[f"{SCORE}_aug"].to_numpy(float))

        for (db, sb), g in m.groupby(["dbin", "sbin"]):
            y = g[LABEL].to_numpy(int)
            ro_b, pr_b = score(y, g[f"{SCORE}_base"].to_numpy(float))
            ro_a, pr_a = score(y, g[f"{SCORE}_aug"].to_numpy(float))
            sub = prot[(prot.dbin == db) & (prot.sbin == sb)]
            rows.append(dict(
                task=task, depth_bin=int(db), size_bin=int(sb),
                n_proteins=len(sub), n_sites=len(y), pos_rate=y.mean(),
                depth_median=sub.depth.median(), size_median=sub["size"].median(),
                d_pct=g.pct_aug.mean() - g.pct_base.mean(),
                d_pct_pos=(g.loc[g[LABEL] == 1, "pct_aug"].mean()
                           - g.loc[g[LABEL] == 1, "pct_base"].mean()),
                d_pct_neg=(g.loc[g[LABEL] == 0, "pct_aug"].mean()
                           - g.loc[g[LABEL] == 0, "pct_base"].mean()),
                d_auroc=ro_a - ro_b, d_auprc=pr_a - pr_b,
                thin=len(y) < args.min_sites))

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False, float_format="%.4f")

    f4 = lambda v: f"{v:+.4f}"
    for task in tasks:
        s = out[out.task == task]
        print(f"\n=== {task}  ({args.family}, {args.train}->{args.eval}) ===")
        piv = s.pivot(index="depth_bin", columns="size_bin", values="d_pct")
        piv["all sizes"] = (s.groupby("depth_bin")
                             .apply(lambda g: np.average(g.d_pct.fillna(0),
                                                         weights=g.n_sites)))
        print("change in mean score percentile, by depth bin (rows) "
              "and candidate-count bin (columns)")
        print(piv.to_string(float_format=f4))
        summ = s.groupby("depth_bin").apply(lambda g: pd.Series({
            "median depth": g.depth_median.iloc[0],
            "pos rate": np.average(g.pos_rate, weights=g.n_sites),
            "d pct (pos)": np.average(g.d_pct_pos.fillna(0), weights=g.n_sites),
            "d pct (neg)": np.average(g.d_pct_neg.fillna(0), weights=g.n_sites),
            "dAUPRC in-bin": np.average(g.d_auprc.fillna(0), weights=g.n_sites),
        }))
        print(summ.to_string(float_format=lambda v: f"{v:+.4f}"))
        thin = int(s.thin.sum())
        if thin:
            print(f"  {thin} of {len(s)} cells hold fewer than {args.min_sites} sites")

    print(f"\nwrote {args.out}")
    print("\nRead the table down a column, not across: a depth effect is a "
          "gradient in percentile that survives within every candidate-count "
          "stratum. A gradient present only in the 'all sizes' column is a size "
          "effect. The shortcut account predicts that the shallow bins gain "
          "percentile while carrying the lower positive rate, and that the "
          "gradient is absent when train=rebuilt.")


if __name__ == "__main__":
    sys.exit(main())
