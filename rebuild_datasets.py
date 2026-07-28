#!/usr/bin/env python3
"""
rebuild_datasets.py — rebuild the PTM datasets changing ONE variable.

DeepMVP's protocol sampled negatives only from proteins carrying at least ten
PTM sites. That rule makes every protein below the threshold contribute
positives and nothing else, so protein identity alone determines the label for a
large share of sites. This script rebuilds the same datasets with the threshold
removed and everything else held fixed:

  * positives           : taken verbatim from retrain/*.tsv (PTMAtlas)
  * candidate negatives : every residue of the target type, in the same proteins,
                          from the same 2019-02-14 Swiss-Prot reference
  * sampling            : uniform at random at a fixed negative:positive ratio,
                          with no per-protein threshold

Run with --stats first to see what the reconstruction looks like before writing
anything.

Usage:
    python rebuild_datasets.py --stats
    python rebuild_datasets.py --ratio 5 --out rebuilt
    python rebuild_datasets.py --ratio 5 --min-sites 10 --out rebuilt_replica
"""

import argparse
import collections
import os
import random
import sys
from itertools import groupby

import pandas as pd

FASTA = "/home/FCAM/juli/HRP/deepmvp/DeepMVP/data/swiss_prot_human_20190214.fasta"
RETRAIN = "/home/FCAM/juli/HRP/retrain"
WINDOW = 31

# PTM -> residues that are candidate sites
TARGETS = {
    "acetylation_k": "K", "glycosylation_n": "N", "methylation_k": "KR",
    "methylation_r": "R", "phosphorylation_st": "ST", "phosphorylation_y": "Y",
    "sumoylation_k": "K", "ubiquitination_k": "K",
}


def read_fasta(path):
    d = {}
    with open(path) as fh:
        for is_hdr, grp in groupby(fh, lambda l: l.startswith(">")):
            if is_hdr:
                acc = next(grp).split("|")[1]
            else:
                d[acc] = "".join(x.strip() for x in grp)
    return d


def window(seq, pos, size=WINDOW):
    """1-based pos; pad with X to a fixed odd length."""
    half = size // 2
    i = pos - 1
    left = seq[max(0, i - half):i].rjust(half, "X")
    right = seq[i + 1:i + 1 + half].ljust(half, "X")
    return left + seq[i] + right


def is_sequon(seq, pos):
    """N-X-S/T with X != P, matching the N!P[ST] filter used for N-glycosylation."""
    i = pos - 1
    if i + 2 >= len(seq):
        return False
    return seq[i + 1] != "P" and seq[i + 2] in "ST"


# DeepMVP trains ONE model per merged PTM type; the >=10-site threshold is
# evaluated on that merged set, not on the individual residue task.
THRESHOLD_GROUP = {
    "acetylation_k": ["acetylation_k"],
    "glycosylation_n": ["glycosylation_n"],
    "methylation_k": ["methylation_k"],          # already the merged K+R set
    "methylation_r": ["methylation_k"],
    "phosphorylation_st": ["phosphorylation_st", "phosphorylation_y"],
    "phosphorylation_y": ["phosphorylation_st", "phosphorylation_y"],
    "sumoylation_k": ["sumoylation_k"],
    "ubiquitination_k": ["ubiquitination_k"],
}


def load_positives(ptm):
    frames = [pd.read_csv(f"{RETRAIN}/{ptm}_{s}.tsv", sep="\t") for s in ("train", "test")]
    df = pd.concat(frames, ignore_index=True)
    return df[df.y == 1][["protein", "aa", "pos"]]


def enumerate_candidates(fa, ptm, proteins, sequon_negatives=False):
    """All residues of the target type in the given proteins.

    N-glycosylation positives in PTMAtlas passed an N!P[ST] sequon filter, but
    the published negatives did not (only ~7% of them are sequons). Default
    behaviour reproduces that asymmetry; --sequon-negatives removes it.
    """
    res = set(TARGETS[ptm])
    out = []
    for p in proteins:
        s = fa.get(p)
        if s is None:
            continue
        for i, c in enumerate(s):
            if c not in res:
                continue
            pos = i + 1
            if ptm == "glycosylation_n" and sequon_negatives and not is_sequon(s, pos):
                continue
            out.append((p, c, pos))
    return out


def build(ptm, fa, ratio, min_sites, seed, sequon_negatives=False):
    pos_df = load_positives(ptm)
    pos_set = {(r.protein, r.pos) for r in pos_df.itertuples()}
    proteins = set(pos_df.protein)

    if min_sites:
        # count sites over the merged PTM type the model is actually trained on
        merged = pd.concat([load_positives(m) for m in THRESHOLD_GROUP[ptm]],
                           ignore_index=True)
        counts = collections.Counter(merged.protein)
        donors = {p for p in set(merged.protein) if counts[p] >= min_sites}
    else:
        donors = proteins

    cands = [c for c in enumerate_candidates(fa, ptm, donors, sequon_negatives)
             if (c[0], c[2]) not in pos_set]

    rng = random.Random(seed)
    want = int(len(pos_set) * ratio)
    negs = cands if len(cands) <= want else rng.sample(cands, want)

    rows = []
    for p, aa, po in pos_df.itertuples(index=False):
        rows.append((p, aa, po, window(fa[p], po), 1))
    for p, aa, po in negs:
        rows.append((p, aa, po, window(fa[p], po), 0))

    df = pd.DataFrame(rows, columns=["protein", "aa", "pos", "x", "y"])

    per = df.groupby("protein")["y"].agg(["size", "sum"])
    pure_pos = int((per["sum"] == per["size"]).sum())
    pure_sites = int(per.loc[per["sum"] == per["size"], "size"].sum())
    stats = dict(
        ptm=ptm, positives=len(pos_set), candidates=len(cands),
        max_ratio=len(cands) / max(len(pos_set), 1), negatives=len(negs),
        proteins=len(per), pure_pos_proteins=pure_pos,
        pure_pos_pct=100 * pure_pos / len(per),
        pure_pos_site_pct=100 * pure_sites / len(df),
        pos_rate=float(df.y.mean()),
    )
    return df, stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", action="store_true", help="report only, write nothing")
    ap.add_argument("--ratio", type=float, default=5.0,
                    help="negatives per positive")
    ap.add_argument("--min-sites", type=int, default=0,
                    help="replicate the original threshold (use 10 to reproduce it)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--sequon-negatives", action="store_true",
                    help="also require the N!P[ST] sequon for N-glyc negatives")
    ap.add_argument("--out", default="rebuilt")
    ap.add_argument("--ptm", default="all")
    a = ap.parse_args()

    if not os.path.exists(FASTA):
        sys.exit(f"missing reference fasta: {FASTA}")
    fa = read_fasta(FASTA)
    print(f"reference: {len(fa)} proteins\n")

    ptms = list(TARGETS) if a.ptm == "all" else [a.ptm]
    hdr = (f"{'PTM':20s} {'pos':>8s} {'cands':>9s} {'max r':>7s} {'negs':>9s} "
           f"{'pos rate':>9s} {'pure+ prot':>11s} {'pure+ sites':>12s}")
    print(hdr)
    print("-" * len(hdr))

    for ptm in ptms:
        df, s = build(ptm, fa, a.ratio, a.min_sites, a.seed,
                      a.sequon_negatives)
        print(f"{s['ptm']:20s} {s['positives']:8d} {s['candidates']:9d} "
              f"{s['max_ratio']:7.1f} {s['negatives']:9d} {s['pos_rate']:9.4f} "
              f"{s['pure_pos_pct']:10.1f}% {s['pure_pos_site_pct']:11.1f}%")
        if not a.stats:
            os.makedirs(a.out, exist_ok=True)
            df.to_csv(f"{a.out}/{ptm}_all.tsv", sep="\t", index=False)

    print("\npure+ prot  = share of proteins whose sites are all positive")
    print("pure+ sites = share of sites living on such proteins (the leak channel)")
    if a.stats:
        print("\n--stats given: nothing written. Re-run without it to emit TSVs.")
    else:
        print(f"\nwrote {len(ptms)} files to {a.out}/ (unsplit; run CD-HIT next)")


if __name__ == "__main__":
    main()
