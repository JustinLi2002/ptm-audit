#!/usr/bin/env python3
"""Independently recompute the 2 x 2 cross-evaluation behind Tables 2, 3, S3
and S4, from the per-site predictions.

This deliberately does NOT call crosseval_summary.py.  Rerunning the original
script would only show it is deterministic; it would not catch an error in the
script itself, which is the failure mode already found twice in this audit
(the AUPRC agreement claim, and the ICC fold-increase).  The point is a second
implementation that can disagree.

Naming convention (verified against 672 files in pdisjoint_runs_v2):

    {task}__{train}__{cond}__split{n}[__{source}]__on_{eval}.pred.tsv

    train  : replica (threshold-sampled) | rebuilt (unrestricted)
    cond   : baseline | ppi (real feature) | shuffled (permuted feature)
    source : absent -> node2vec interaction embedding
             esm    -> ESM-2
             prott5 -> ProtT5-XL-U50
    eval   : replica (threshold-drawn negatives) | rebuilt (natural negatives)

Delta is computed per (task, train, eval, split) as
    metric(feature model) - metric(sequence-only baseline)
and then averaged over splits, so the baseline is always the one from the
matching cell.  Columns are asserted, never inferred by name.

Usage
-----
    python crosseval_verify.py --root ~/HRP/pdisjoint_runs_v2 --dry-run
    python crosseval_verify.py --root ~/HRP/pdisjoint_runs_v2 -o crosseval_verify.tsv
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

# The four cells of the 2 x 2, in the order the manuscript presents them.
CELLS = [("replica", "replica", "threshold train, threshold eval"),
         ("replica", "rebuilt", "threshold train, natural eval"),
         ("rebuilt", "replica", "unrestricted train, threshold eval"),
         ("rebuilt", "rebuilt", "unrestricted train, natural eval")]


def family(cond, source):
    if cond == "baseline":
        return "baseline"
    src = {None: "interaction", "esm": "esm2", "prott5": "prott5"}[source]
    return src if cond == "ppi" else src + "_perm"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.expanduser("~/HRP/pdisjoint_runs_v2"))
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-o", "--out", default="crosseval_verify.tsv")
    args = ap.parse_args()

    parsed, unmatched = [], []
    for p in sorted(glob.glob(os.path.join(args.root, "*.pred.tsv"))):
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
        print(f"matched {len(parsed)} / {len(parsed) + len(unmatched)}")
        print(df.groupby(["family", "train", "eval"]).size().to_string())
        for u in unmatched[:10]:
            print("  UNMATCHED:", u)
        return
    if unmatched:
        raise SystemExit(f"{len(unmatched)} files did not parse; use --dry-run")

    rows = []
    for d in parsed:
        df = pd.read_csv(d["path"], sep="\t")
        missing = {PROTEIN, SITE, LABEL, SCORE} - set(df.columns)
        assert not missing, f"{d['path']}: missing {missing}"
        y = df[LABEL].to_numpy(int)
        s = df[SCORE].to_numpy(float)
        assert set(np.unique(y)) <= {0, 1}, f"{d['path']}: labels not binary"
        rows.append({k: d[k] for k in ("task", "train", "eval", "split", "family")}
                    | dict(auroc=roc_auc_score(y, s),
                           auprc=average_precision_score(y, s),
                           n_sites=len(y), pos_rate=y.mean()))

    long = pd.DataFrame(rows)
    key = ["task", "train", "eval", "split"]

    base = (long[long.family == "baseline"].set_index(key)
            [["auroc", "auprc", "n_sites", "pos_rate"]]
            .rename(columns=lambda c: c + "_base"))
    out = (long[long.family != "baseline"].set_index(key)
           .join(base, how="inner").reset_index())
    assert len(out) == len(long[long.family != "baseline"]), \
        "some feature runs had no matching baseline in the same cell"
    # sanity: the baseline and the feature model must have scored the same sites
    bad = out[out.n_sites != out.n_sites_base]
    assert bad.empty, f"site counts differ from baseline in {len(bad)} rows"

    out["d_auroc"] = out.auroc - out.auroc_base
    out["d_auprc"] = out.auprc - out.auprc_base
    out.to_csv(args.out, sep="\t", index=False, float_format="%.4f")

    agg = (out.groupby(["family", "train", "eval", "task"])
              [["d_auroc", "d_auprc", "auroc_base", "auroc"]].mean().reset_index())
    agg["_o"] = agg.task.map({t: i for i, t in enumerate(ORDER)})

    for fam in ["interaction", "esm2", "prott5",
                "interaction_perm", "esm2_perm", "prott5_perm"]:
        sub = agg[agg.family == fam]
        if sub.empty:
            continue
        print(f"\n{'='*78}\n{fam}\n{'='*78}")
        piv = sub.pivot_table(index="task", columns=["train", "eval"],
                              values="d_auroc")
        piv = piv.reindex([t for t in ORDER if t in piv.index])
        piv = piv[[(tr, ev) for tr, ev, _ in CELLS if (tr, ev) in piv.columns]]
        piv.columns = [lab for tr, ev, lab in CELLS if (tr, ev) in
                       set(zip(sub.train, sub['eval']))]
        print("dAUROC")
        print(piv.to_string(float_format=lambda v: f"{v:+.4f}"))
        print("mean: " + "  ".join(f"{c}={piv[c].mean():+.4f}" for c in piv.columns))
        print("negative in: " + "  ".join(
            f"{c}={int((piv[c] < 0).sum())}/{len(piv)}" for c in piv.columns))

        pivp = sub.pivot_table(index="task", columns=["train", "eval"],
                               values="d_auprc")
        pivp = pivp.reindex(piv.index)
        pivp = pivp[[(tr, ev) for tr, ev, _ in CELLS if (tr, ev) in pivp.columns]]
        pivp.columns = piv.columns
        print("\ndAUPRC")
        print(pivp.to_string(float_format=lambda v: f"{v:+.4f}"))
        print("mean: " + "  ".join(f"{c}={pivp[c].mean():+.4f}" for c in pivp.columns))

    print(f"\n{'='*78}\nreal minus permuted, dAUROC (same cell, same split)\n{'='*78}")
    for fam in ["interaction", "esm2", "prott5"]:
        r = agg[agg.family == fam].set_index("task").d_auroc
        p = agg[agg.family == fam + "_perm"].set_index("task").d_auroc
        if r.empty or p.empty:
            continue
        # recompute per cell rather than collapsing
        rr = (out[out.family == fam].groupby(["train", "eval", "task"]).d_auroc.mean())
        pp = (out[out.family == fam + "_perm"]
              .groupby(["train", "eval", "task"]).d_auroc.mean())
        diff = (rr - pp).unstack(["train", "eval"])
        diff = diff.reindex([t for t in ORDER if t in diff.index])
        print(f"\n{fam}")
        print(diff.to_string(float_format=lambda v: f"{v:+.4f}"))
        print("real ahead in: " + "  ".join(
            f"{c}={int((diff[c] > 0).sum())}/{len(diff)}" for c in diff.columns))

    print(f"\nwrote {args.out}  ({len(out)} rows)")
    print("\nCompare against the manuscript:")
    print("  Table 3 (language model): threshold/threshold +0.0437,")
    print("      threshold/natural -0.1008 AUROC and -0.133 AUPRC,")
    print("      unrestricted/natural +0.0376; real > permuted in 2/8 tasks")
    print("  Table 2 (interaction):    threshold/threshold -0.0280,")
    print("      threshold/natural -0.1250, unrestricted +0.0200")


if __name__ == "__main__":
    sys.exit(main())
