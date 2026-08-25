#!/usr/bin/env python3
"""Share of deployment negatives lying on proteins below the donor threshold.

w_k is stratum k's share of the negatives a predictor actually meets, and needs
only a candidate-residue count per protein and that protein's annotation depth
-- no labels, no negatives from the benchmark, no retraining. The quantity
reported here is sum_{d<T} w_k: the part of the deployment negative population
about which a threshold-sampled benchmark carries no information at all.

Depth is the merged-type count the donor rule is actually evaluated on, taken
from the unrestricted reconstruction so it does not depend on the rule.

Two universes are reported, because they answer different questions:

  seen   proteins carrying at least one site of the task -- the population the
         reconstruction covers, and a conservative floor
  all    every protein in the reference proteome, which is the population a
         predictor is applied to: most proteins it meets carry no annotated
         site at all, so they sit below any positive threshold by definition
"""
import os
from collections import defaultdict
from itertools import groupby

import pandas as pd

BASE = os.environ.get("PTM_AUDIT_BASE", "/home/FCAM/juli/HRP")
FASTA = f"{BASE}/deepmvp/DeepMVP/data/swiss_prot_human_20190214.fasta"
REBUILT = f"{BASE}/rebuilt"
T = 10

TARGETS = {"acetylation_k": "K", "glycosylation_n": "N", "methylation_k": "KR",
           "methylation_r": "R", "phosphorylation_st": "ST",
           "phosphorylation_y": "Y", "sumoylation_k": "K",
           "ubiquitination_k": "K"}
MERGED = {"acetylation_k": ["acetylation_k"],
          "glycosylation_n": ["glycosylation_n"],
          "methylation_k": ["methylation_k"],
          "methylation_r": ["methylation_k"],
          "phosphorylation_st": ["phosphorylation_st", "phosphorylation_y"],
          "phosphorylation_y": ["phosphorylation_st", "phosphorylation_y"],
          "sumoylation_k": ["sumoylation_k"],
          "ubiquitination_k": ["ubiquitination_k"]}
ORDER = ["phosphorylation_st", "phosphorylation_y", "acetylation_k",
         "methylation_k", "methylation_r", "sumoylation_k",
         "ubiquitination_k", "glycosylation_n"]
SHORT = {"phosphorylation_st": "Phosphorylation S/T", "phosphorylation_y": "Phosphorylation Y",
         "acetylation_k": "Acetylation K", "methylation_k": "Methylation K/R",
         "methylation_r": "Methylation R", "sumoylation_k": "Sumoylation K",
         "ubiquitination_k": "Ubiquitination K", "glycosylation_n": "N-Glycosylation N"}


def read_fasta(path):
    d = {}
    with open(path) as fh:
        for is_hdr, grp in groupby(fh, lambda l: l.startswith(">")):
            if is_hdr:
                acc = next(grp).split("|")[1]
            else:
                d[acc] = "".join(x.strip() for x in grp)
    return d


def positives():
    """{task: {protein: set(pos)}} from the unrestricted reconstruction."""
    out = {}
    for t in TARGETS:
        df = pd.read_csv(f"{REBUILT}/{t}_all.tsv", sep="\t", usecols=["protein", "pos", "y"])
        df = df[df.y == 1]
        d = defaultdict(set)
        for p, pos in zip(df.protein, df.pos):
            d[p].add(int(pos))
        out[t] = d
    return out


def main():
    fa = read_fasta(FASTA)
    pos = positives()
    print(f"reference proteome: {len(fa)} sequences;  threshold T = {T}\n")
    hdr = (f"{'task':22s} {'seen: unobs w':>14s} {'n prot':>8s} "
           f"{'all: unobs w':>13s} {'n prot':>8s}")
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for t in ORDER:
        res = set(TARGETS[t])
        depth = defaultdict(int)
        for m in MERGED[t]:
            for p, s in pos[m].items():
                depth[p] += len(s)
        seen_prot = set(pos[t])
        num_seen = den_seen = num_all = den_all = 0
        for p, seq in fa.items():
            cand = sum(1 for c in seq if c in res) - len(pos[t].get(p, ()))
            if cand <= 0:
                continue
            below = depth.get(p, 0) < T
            den_all += cand
            if below:
                num_all += cand
            if p in seen_prot:
                den_seen += cand
                if below:
                    num_seen += cand
        ws = num_seen / den_seen if den_seen else float("nan")
        wa = num_all / den_all if den_all else float("nan")
        rows.append((t, ws, wa))
        print(f"{SHORT[t]:22s} {ws:14.3f} {len(seen_prot):8d} "
              f"{wa:13.3f} {sum(1 for p in fa if sum(1 for c in fa[p] if c in res) > 0):8d}")
    print()
    print(f"  seen universe : {min(r[1] for r in rows):.3f} - {max(r[1] for r in rows):.3f}")
    print(f"  all  universe : {min(r[2] for r in rows):.3f} - {max(r[2] for r in rows):.3f}")


if __name__ == "__main__":
    main()
