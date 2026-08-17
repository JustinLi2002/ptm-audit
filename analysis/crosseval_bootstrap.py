#!/usr/bin/env python3
"""Paired cluster bootstrap for the cross-evaluation cells of Tables 2, 3 and S11.

For each (feature family, training construction, evaluation construction, task)
cell, resamples whole test proteins with replacement and recomputes the
baseline and the augmented metric on the same resampled proteins within each
iteration, then takes the percentile interval of their difference. Pairing
matters: the two models see an identical protein sample in every iteration, so
the interval on the difference is much narrower than the difference of two
marginal intervals would suggest. This is the procedure already used for the
identity baseline, applied to the cells that carry the headline result.

The three partitions hold different test proteins, so they cannot be pooled
into one urn. Within each iteration every partition is resampled separately and
the resampled partitions are then combined into a single metric, so the
resulting interval carries both protein-sampling and partition variance. That
is wider than a within-partition interval and is the honest quantity here,
since the reported point estimates are means over partitions.

Naming convention (verified against 672 files):

    {task}__{train}__{cond}__split{n}[__{source}]__on_{eval}.pred.tsv

Usage
-----
    python crosseval_bootstrap.py --root ~/HRP/pdisjoint_runs_v2 --dry-run
    python crosseval_bootstrap.py --root ~/HRP/pdisjoint_runs_v2 \
        -B 2000 -o crosseval_ci.tsv
"""
import argparse
import glob
import hashlib
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

FAMILIES = ["interaction", "esm2", "prott5",
            "interaction_perm", "esm2_perm", "prott5_perm"]


def family(cond, source):
    if cond == "baseline":
        return "baseline"
    src = {None: "interaction", "esm": "esm2", "prott5": "prott5"}[source]
    return src if cond == "ppi" else src + "_perm"


def load(path):
    df = pd.read_csv(path, sep="\t")
    missing = {PROTEIN, SITE, LABEL, SCORE} - set(df.columns)
    assert not missing, f"{path}: missing {missing}"
    return df


def metrics(y, s):
    return roc_auc_score(y, s), average_precision_score(y, s)


def bootstrap_cell(parts, B, rng):
    """parts: list of (protein_codes, y, s_base, s_aug), one per partition.

    Returns the point estimates and percentile intervals of the paired
    differences. Within an iteration each partition is resampled over its own
    proteins; the resampled partitions are then concatenated and scored once.
    """
    # index sites by protein, per partition, once
    prepared = []
    for codes, y, sb, sa in parts:
        order = np.argsort(codes, kind="stable")
        codes, y, sb, sa = codes[order], y[order], sb[order], sa[order]
        uniq, start = np.unique(codes, return_index=True)
        bounds = np.append(start, len(codes))
        # sizes and offsets let a resample be built with one repeat + one
        # arange instead of a Python loop over proteins
        sizes = np.diff(bounds)
        prepared.append((bounds[:-1], sizes, len(uniq), y, sb, sa))

    d_auroc = np.empty(B)
    d_auprc = np.empty(B)
    kept = 0
    for b in range(B):
        idx = []
        for offs, sizes, n_prot, y, sb, sa in prepared:
            pick = rng.integers(0, n_prot, n_prot)
            take = sizes[pick]
            starts = np.repeat(offs[pick], take)
            within = np.arange(take.sum()) - np.repeat(np.cumsum(take) - take, take)
            idx.append(starts + within)
        ys = np.concatenate([p[3][i] for p, i in zip(prepared, idx)])
        if ys.min() == ys.max():          # degenerate resample, no both classes
            continue
        bs = np.concatenate([p[4][i] for p, i in zip(prepared, idx)])
        as_ = np.concatenate([p[5][i] for p, i in zip(prepared, idx)])
        ro_b, pr_b = metrics(ys, bs)
        ro_a, pr_a = metrics(ys, as_)
        d_auroc[kept] = ro_a - ro_b
        d_auprc[kept] = pr_a - pr_b
        kept += 1

    y_all = np.concatenate([p[3] for p in prepared])
    b_all = np.concatenate([p[4] for p in prepared])
    a_all = np.concatenate([p[5] for p in prepared])
    ro_b, pr_b = metrics(y_all, b_all)
    ro_a, pr_a = metrics(y_all, a_all)

    q = lambda v: np.percentile(v[:kept], [2.5, 97.5])
    lo_ro, hi_ro = q(d_auroc)
    lo_pr, hi_pr = q(d_auprc)
    return dict(d_auroc=ro_a - ro_b, auroc_lo=lo_ro, auroc_hi=hi_ro,
                d_auprc=pr_a - pr_b, auprc_lo=lo_pr, auprc_hi=hi_pr,
                n_boot=kept, n_sites=len(y_all),
                n_proteins=sum(p[2] for p in prepared))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.expanduser("~/HRP/pdisjoint_runs_v2"))
    ap.add_argument("-B", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--families", nargs="*", default=FAMILIES)
    ap.add_argument("--shard", type=int, default=None,
                    help="0-based index of a single cell to run (for job arrays)")
    ap.add_argument("--n-shards", type=int, default=None,
                    help="if set with --shard, run every n-th cell instead of one")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-o", "--out", default="crosseval_ci.tsv")
    args = ap.parse_args()

    files = {}
    for p in sorted(glob.glob(os.path.join(args.root, "*.pred.tsv"))):
        m = re.match(PATTERN, os.path.basename(p))
        if not m:
            raise SystemExit(f"did not parse: {os.path.basename(p)}")
        d = m.groupdict()
        key = (d["task"], d["train"], d["eval"], family(d["cond"], d["source"]))
        files.setdefault(key, {})[d["split"]] = p

    cells = [(t, tr, ev, f)
             for f in args.families
             for tr in ("replica", "rebuilt")
             for ev in ("replica", "rebuilt")
             for t in ORDER
             if (t, tr, ev, f) in files]

    if args.dry_run:
        print(f"{len(cells)} cells, {args.B} resamples each")
        for c in cells[:5]:
            print("  ", c, sorted(files[c]))
        print("   ...")
        return

    if args.shard is not None:
        if args.n_shards:
            cells = cells[args.shard::args.n_shards]
        else:
            cells = [cells[args.shard]]

    rows = []
    for n, (task, tr, ev, fam) in enumerate(cells, 1):
        # seed from the cell identity, not from position, so a shard gives the
        # same numbers as the corresponding rows of a single whole-set run
        # hashlib, not hash(): PYTHONHASHSEED randomises str hashing per
        # process, so hash() would give a different stream in every run and a
        # shard would not reproduce the corresponding row of a whole-set run
        tag = f"{args.seed}|{task}|{tr}|{ev}|{fam}".encode()
        rng = np.random.default_rng(
            int.from_bytes(hashlib.sha256(tag).digest()[:8], "little"))
        splits = files[(task, tr, ev, fam)]
        base = files[(task, tr, ev, "baseline")]
        assert set(splits) == set(base), f"{task}/{tr}/{ev}/{fam}: split mismatch"
        parts = []
        for s in sorted(splits):
            a, b = load(base[s]), load(splits[s])
            key = [PROTEIN, SITE]
            m = a.merge(b, on=key + [LABEL], suffixes=("_base", "_aug"),
                        validate="one_to_one")
            assert len(m) == len(a), f"{task}/{tr}/{ev}/{fam}/{s}: merge lost rows"
            parts.append((m[PROTEIN].to_numpy(),
                          m[LABEL].to_numpy(int),
                          m[f"{SCORE}_base"].to_numpy(float),
                          m[f"{SCORE}_aug"].to_numpy(float)))
        r = bootstrap_cell(parts, args.B, rng)
        rows.append(dict(task=task, train=tr, eval=ev, family=fam, **r))
        print(f"  [{n}/{len(cells)}] {fam:17s} {tr}->{ev} {task:20s} "
              f"dAUROC {r['d_auroc']:+.4f} ({r['auroc_lo']:+.4f}, {r['auroc_hi']:+.4f})",
              flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False, float_format="%.4f")
    print(f"\nwrote {args.out}")

    print("\ncells whose 95% interval on dAUROC contains zero:")
    z = out[(out.auroc_lo < 0) & (out.auroc_hi > 0)]
    if z.empty:
        print("  none")
    else:
        print(z[["family", "train", "eval", "task", "d_auroc",
                 "auroc_lo", "auroc_hi"]]
              .to_string(index=False, float_format=lambda v: f"{v:+.4f}"))


if __name__ == "__main__":
    sys.exit(main())
