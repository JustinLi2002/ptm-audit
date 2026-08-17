#!/usr/bin/env python3
"""Merge the per-cell bootstrap shards and check that none is missing or thin.

    python merge_ci_shards.py ~/HRP/ci_shards -o crosseval_ci.tsv
"""
import argparse
import glob
import os
import sys

import pandas as pd

EXPECTED = 192


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("shards", nargs="?", default=os.path.expanduser("~/HRP/ci_shards"))
    ap.add_argument("-o", "--out", default="crosseval_ci.tsv")
    ap.add_argument("--min-boot", type=int, default=1900,
                    help="flag cells where too many resamples were degenerate")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.shards, "cell_*.tsv")))
    if len(files) != EXPECTED:
        have = {int(os.path.basename(f)[5:8]) for f in files}
        missing = sorted(set(range(EXPECTED)) - have)
        print(f"{len(files)}/{EXPECTED} shards present.")
        print("missing:", missing)
        print("resubmit with:  sbatch --array=" +
              ",".join(map(str, missing)) + " crosseval_ci.slurm")
        if not files:
            return 1

    out = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    out.to_csv(args.out, sep="\t", index=False, float_format="%.4f")
    print(f"wrote {args.out}  ({len(out)} cells)")

    thin = out[out.n_boot < args.min_boot]
    if not thin.empty:
        print(f"\ncells with fewer than {args.min_boot} usable resamples "
              "(degenerate draws were skipped; treat their intervals with care):")
        print(thin[["family", "train", "eval", "task", "n_boot", "n_proteins"]]
              .to_string(index=False))

    print("\n95% interval on dAUROC contains zero:")
    z = out[(out.auroc_lo < 0) & (out.auroc_hi > 0)]
    print("  none" if z.empty else
          z[["family", "train", "eval", "task", "d_auroc", "auroc_lo", "auroc_hi"]]
          .to_string(index=False, float_format=lambda v: f"{v:+.4f}"))

    print("\n95% interval on dAUPRC contains zero:")
    z = out[(out.auprc_lo < 0) & (out.auprc_hi > 0)]
    print("  none" if z.empty else
          z[["family", "train", "eval", "task", "d_auprc", "auprc_lo", "auprc_hi"]]
          .to_string(index=False, float_format=lambda v: f"{v:+.4f}"))

    print("\nper-cell mean over tasks, dAUROC:")
    g = (out.groupby(["family", "train", "eval"])
            [["d_auroc", "d_auprc"]].mean())
    print(g.to_string(float_format=lambda v: f"{v:+.4f}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
