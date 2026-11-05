#!/usr/bin/env python3
"""Ratio-matched training construction for the mouse replication (PREREG_mouse.md, "Constructions").

The rule of ratio_matched_v4/build_matched.py, unchanged. For each split seed, and for each of the train and
validation partitions, it takes the unrestricted positives plus a uniform subsample of the unrestricted
negatives, of the threshold-sampled negative count. Test rows are not written.

Also reports the pure-positive share of training sites in T, U and matched.
"""
import hashlib
import os
import sys
from pathlib import Path

import pandas as pd

BASE = Path(os.environ["PTM_AUDIT_BASE"])
assert BASE.name == "mouse_v1", BASE
OUT = BASE / "matched"
TASKS = ["phosphorylation_st"]
SPLITS = [0, 1, 2]
PARTS = {"train": 1000, "val": 2000}        # sampling seed = offset + split seed


def pure_pos_share(d):
    """Share of sites on proteins whose retained sites are all positive."""
    return float((d.groupby("protein").y.transform("min") == 1).mean())


def sites(d):
    return set(zip(d.protein, d.pos))


rows, files = [], []
for task in TASKS:
    T = pd.read_csv(BASE / "replica" / f"{task}_all.tsv", sep="\t")
    U = pd.read_csv(BASE / "rebuilt" / f"{task}_all.tsv", sep="\t")
    for s in SPLITS:
        sp = pd.read_csv(BASE / "pdisjoint" / f"split_seed{s}.csv")
        sp = dict(zip(sp.accession, sp.split))
        out = []
        for part, off in PARTS.items():
            t = T[T.protein.map(sp) == part]
            u = U[U.protein.map(sp) == part]
            pt, pu = sites(t[t.y == 1]), sites(u[u.y == 1])
            neg_t = int((t.y == 0).sum())
            pool = u[u.y == 0]
            n_neg = neg_t if pt == pu else int(round(neg_t * len(pu) / len(pt)))
            if n_neg > len(pool):
                sys.exit(f"{task} split{s} {part}: needs {n_neg} negatives, unrestricted pool has {len(pool)}")
            m = pd.concat([u[u.y == 1], pool.sample(n=n_neg, random_state=off + s)]).sort_index()
            out.append(m)
            rows.append(dict(task=task, split=s, part=part, pos_T=len(pt), pos_U=len(pu),
                             pos_shared=len(pt & pu), neg_T=neg_t, neg_U=len(pool), neg_matched=n_neg,
                             rate_T=t.y.mean(), rate_U=u.y.mean(), rate_matched=m.y.mean(),
                             prot_T=t.protein.nunique(), prot_U=u.protein.nunique(),
                             prot_matched=m.protein.nunique(), purepos_T=pure_pos_share(t),
                             purepos_U=pure_pos_share(u), purepos_matched=pure_pos_share(m)))
        d = OUT / f"split{s}"
        d.mkdir(parents=True, exist_ok=True)
        f = d / f"{task}_all.tsv"
        pd.concat(out).to_csv(f, sep="\t", index=False)
        files.append(f)

summ = pd.DataFrame(rows)
summ.to_csv(OUT / "matched_summary.tsv", sep="\t", index=False)
with open(OUT / "manifest.sha256", "w") as fh:
    for f in files:
        fh.write(f"{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(OUT)}\n")
pd.set_option("display.width", 250)
print(summ.round(4).to_string())
