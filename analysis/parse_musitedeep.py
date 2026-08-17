#!/usr/bin/env python3
"""Convert a MusiteDeep-style annotated FASTA into (protein, pos, y) rows.

MusiteDeep releases full-length sequences with modified residues marked by a
character inserted after them, by default '#':

    >sp|Q4R5L1|AATC_MACFA Aspartate aminotransferase, cytoplasmic ...
    MAPPSVFSEVPQ...VFS#AAGFKD...

Every residue of the target type is a candidate. Marked ones are positives and
the rest are negatives, so a protein contributes negatives whatever its
annotation depth. That is the unrestricted construction, and reading these
files is how one checks whether a benchmark applies a donor threshold at all
rather than taking the method section's word for it.

The residue alphabet per modification is taken from the directory name where
possible and can be overridden.

    python parse_musitedeep.py --fasta <file> --residues ST -o phospho_st.tsv
    python parse_musitedeep.py --root MusiteDeep_web/MusiteDeep/testdata \
        --only train -o musitedeep_all.tsv
"""
import argparse
import glob
import os
import re
import sys

import pandas as pd

# directory name -> residues that are candidates for that modification
RESIDUES = {
    "ST": "ST", "Y": "Y", "K": "K", "N": "N", "R": "R", "C": "C", "Q": "Q",
    "Phosphorylation": "STY", "N-linked": "N", "O-linked": "ST",
    "Ubiquitination": "K", "SUMOylation": "K", "Acetylation": "K",
    "Methylation": "KR", "Pyrrolidone": "Q", "Palmitoylation": "C",
    "Hydroxylation": "PK", "Glutathionylation": "C", "Nitrosylation": "C",
}
ACC = re.compile(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b")


def residues_for(path, override):
    if override:
        return override.upper()
    parts = os.path.normpath(path).split(os.sep)
    for p in reversed(parts):
        if p in RESIDUES:
            return RESIDUES[p]
    raise SystemExit(f"cannot infer candidate residues from {path}; pass --residues")


def parse(path, residues, mark="#"):
    rows, name, seq = [], None, []

    def flush():
        if name is None:
            return
        s = "".join(seq)
        pos = 0                      # index into the unmarked sequence, 1-based
        i = 0
        while i < len(s):
            ch = s[i]
            if ch == mark:           # a mark refers to the residue before it
                i += 1
                continue
            pos += 1
            if ch in residues:
                marked = i + 1 < len(s) and s[i + 1] == mark
                rows.append((name, pos, int(marked)))
            i += 1

    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            flush()
            m = ACC.search(line)
            name = m.group(1) if m else line[1:].split()[0]
            seq = []
        elif line:
            seq.append(line)
    flush()
    if not rows:
        raise SystemExit(f"{path}: no candidate residues found for {residues!r}")
    return pd.DataFrame(rows, columns=["protein", "pos", "y"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta")
    ap.add_argument("--root")
    ap.add_argument("--only", default=None, help="substring the filename must contain")
    ap.add_argument("--residues", default=None)
    ap.add_argument("--mark", default="#")
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()

    if args.fasta:
        files = [args.fasta]
    elif args.root:
        files = [p for p in glob.glob(os.path.join(args.root, "**", "*.fasta"),
                                      recursive=True)
                 if "annotated" in os.path.basename(p)]
        if args.only:
            files = [p for p in files if args.only in os.path.basename(p)]
    else:
        raise SystemExit("give --fasta or --root")
    if not files:
        raise SystemExit("no annotated fasta matched")

    frames = []
    for p in sorted(files):
        res = residues_for(p, args.residues)
        df = parse(p, res, args.mark)
        tag = os.path.relpath(p, args.root or os.path.dirname(p))
        df["file"] = tag
        frames.append(df)
        print(f"  {tag:58s} residues={res:4s} "
              f"{df.protein.nunique():6d} proteins {len(df):9d} sites "
              f"pos_rate={df.y.mean():.4f}")

    out = pd.concat(frames, ignore_index=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"\nwrote {args.out}  ({len(out):,} rows)")
    print("Feed it to benchmark_homogeneity.py with --data on this file.")


if __name__ == "__main__":
    sys.exit(main())
