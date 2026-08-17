#!/usr/bin/env python3
"""Label-homogeneity statistics for any PTM site benchmark.

The manuscript states that the donor-threshold convention is not specific to
one resource. That claim currently rests on reading other papers' methods
rather than on measuring their data. This computes, for any dataset that can be
reduced to (protein, position, label) rows, the two quantities the argument
turns on:

    homogeneous proteins   proteins whose sites all carry one label
    pure-positive share    the share of sites that sit on all-positive proteins

Both are reported over the training partition where one is distinguishable, and
over the whole dataset otherwise. Nothing here needs a trained model, so a
dataset can be characterised from its released files alone.

Input formats
-------------
The script tries, in order: a tsv/csv with recognisable column names; a FASTA
of windows whose headers carry an accession and a position; and a directory of
per-task files matching *_train.tsv / *_test.tsv or *_all.tsv.

    python benchmark_homogeneity.py --data ~/HRP/retrain --name "DeepMVP (earlier)"
    python benchmark_homogeneity.py --data ~/musitedeep/data --name MusiteDeep
    python benchmark_homogeneity.py --data ... --name ... -o homogeneity.tsv

Run it on our own reconstructions first: the numbers must reproduce Table 1,
which is what makes the figures for the other benchmarks comparable.
"""
import argparse
import glob
import os
import re
import sys

import pandas as pd

PROTEIN_COLS = ["protein", "prot", "uniprot", "accession", "acc", "protein_id",
                "entry", "id"]
POS_COLS = ["pos", "position", "site", "residue", "idx", "index"]
LABEL_COLS = ["y", "label", "target", "class", "is_positive", "positive"]


def pick(cols, cands, what, path):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c in low:
            return low[c]
    raise SystemExit(f"{path}: no {what} column; saw {list(cols)}")


def from_table(path):
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    p = pick(df.columns, PROTEIN_COLS, "protein", path)
    y = pick(df.columns, LABEL_COLS, "label", path)
    try:
        s = pick(df.columns, POS_COLS, "position", path)
    except SystemExit:
        s = None
        df = df.assign(_pos=range(len(df)))
    out = df[[p, s or "_pos", y]].copy()
    out.columns = ["protein", "pos", "y"]
    out["y"] = out.y.astype(int)
    return out


def from_fasta(path):
    """Window FASTA whose header carries an accession and a position, with the
    label encoded as a 1/0 or pos/neg token somewhere in the header."""
    rows = []
    acc = re.compile(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])")
    for line in open(path):
        if not line.startswith(">"):
            continue
        h = line[1:].strip()
        m = acc.search(h)
        nums = re.findall(r"\b\d+\b", h)
        lab = 1 if re.search(r"\b(1|pos|positive)\b", h, re.I) else 0
        if m and nums:
            rows.append((m.group(1), int(nums[-1]), lab))
    if not rows:
        raise SystemExit(f"{path}: no headers parsed; inspect the format")
    return pd.DataFrame(rows, columns=["protein", "pos", "y"])


def load(path):
    if os.path.isfile(path):
        rd = from_fasta if path.endswith((".fa", ".fasta")) else from_table
        return {os.path.basename(path): rd(path)}
    found = {}
    for p in sorted(glob.glob(os.path.join(path, "*"))):
        if not p.endswith((".tsv", ".csv", ".txt", ".fa", ".fasta")):
            continue
        name = os.path.basename(p)
        try:
            found[name] = from_fasta(p) if p.endswith((".fa", ".fasta")) else from_table(p)
        except SystemExit as e:
            print(f"  skipped {name}: {e}", file=sys.stderr)
    if not found:
        raise SystemExit(f"nothing loadable under {path}")
    return found


def stats(df):
    g = df.groupby("protein").y
    n, mean = g.size(), g.mean()
    homo = (mean == 0) | (mean == 1)
    pure = mean == 1
    sites = n.sum()
    return dict(
        n_proteins=len(n), n_sites=int(sites),
        pos_rate=round(float(df.y.mean()), 4),
        homogeneous_proteins_pct=round(100 * float(homo.mean()), 1),
        homogeneous_sites_pct=round(100 * float(n[homo].sum() / sites), 1),
        pure_positive_sites_pct=round(100 * float(n[pure].sum() / sites), 1),
        median_sites_per_protein=int(n.median()),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--name", required=True, help="benchmark name for the table")
    ap.add_argument("--only", default=None,
                    help="substring a filename must contain, e.g. 'train'")
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

    frames = load(args.data)
    if args.only:
        frames = {k: v for k, v in frames.items() if args.only in k}
        if not frames:
            raise SystemExit(f"no file contains {args.only!r}")

    rows = []
    for name, df in sorted(frames.items()):
        before = len(df)
        df = df.drop_duplicates(["protein", "pos"])
        if len(df) != before:
            print(f"  {name}: dropped {before - len(df)} duplicate (protein, pos) rows",
                  file=sys.stderr)
        rows.append(dict(benchmark=args.name, file=name, **stats(df)))

    out = pd.DataFrame(rows)
    cols = ["benchmark", "file", "n_proteins", "n_sites", "pos_rate",
            "homogeneous_proteins_pct", "homogeneous_sites_pct",
            "pure_positive_sites_pct", "median_sites_per_protein"]
    out = out[cols]
    print(out.to_string(index=False))
    if args.out:
        hdr = not os.path.exists(args.out)
        out.to_csv(args.out, sep="\t", index=False, mode="a", header=hdr)
        print(f"\nappended to {args.out}")

    print(f"\nacross files: homogeneous sites "
          f"{out.homogeneous_sites_pct.min():.1f}-{out.homogeneous_sites_pct.max():.1f}%, "
          f"pure-positive sites "
          f"{out.pure_positive_sites_pct.min():.1f}-{out.pure_positive_sites_pct.max():.1f}%")


if __name__ == "__main__":
    sys.exit(main())
