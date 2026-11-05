#!/usr/bin/env python3
"""Paired protein bootstrap for the cross-evaluation cells of Table 1, matched to the reported statistic.

The point estimate in Table 1 is the mean over the three protein-disjoint splits of the per-split ΔAUROC
(augmented minus sequence-only, each scored within its own split). analysis/crosseval_bootstrap.py instead pooled
the three resampled splits and scored them once, which is a different statistic: for methylation K/R, ProtT5, TT
the split mean is +0.0716 and the pooled value +0.0372. Here every resample computes the reported statistic:

  1. whole proteins are drawn with replacement from the union of the test proteins of the three splits, so a
     protein that lies in the test partition of more than one split carries the same multiplicity in each;
  2. within each split, every site is weighted by its protein's multiplicity, and the sequence-only and augmented
     models are scored on the same weighted sites (weighted Mann-Whitney AUROC, ties counted one half);
  3. the per-split differences are averaged over splits.

The interval is the 2.5th to 97.5th percentile of that mean over resamples. One stream of protein draws per
(task, training construction, test construction) is shared by the three feature families, seeded from the cell
identity so that any subset of cells reproduces the corresponding rows of a full run. The unweighted point
estimate is checked against sklearn on the same files.

Usage: python crosseval_bootstrap_v2.py PRED_DIR OUT_TSV [-B 2000] [--seed 0] [--tasks ...]"""
import argparse
import hashlib
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

TASKS = ["phosphorylation_st", "phosphorylation_y", "acetylation_k", "methylation_k", "methylation_r",
         "sumoylation_k", "ubiquitination_k", "glycosylation_n"]
CELLS = [("replica", "replica"), ("replica", "rebuilt"), ("rebuilt", "replica"), ("rebuilt", "rebuilt")]
FAMILIES = {"interaction": "", "esm2": "__esm", "prott5": "__prott5"}  # file suffix of the augmented model
SPLITS = (0, 1, 2)


class Scored:
    """One model's scores on one split, prepared for repeated weighted AUROC."""

    def __init__(self, y, s):
        order = np.argsort(s, kind="mergesort")
        ss = s[order]
        g_sorted = np.cumsum(np.r_[True, ss[1:] != ss[:-1]]) - 1
        self.g = np.empty(len(s), np.int64)
        self.g[order] = g_sorted
        self.n_groups = int(g_sorted[-1]) + 1
        self.pos = y.astype(bool)

    def auroc(self, w):
        wp = np.bincount(self.g, weights=w * self.pos, minlength=self.n_groups)
        wn = np.bincount(self.g, weights=w * ~self.pos, minlength=self.n_groups)
        below = np.cumsum(wn) - wn
        return float((wp * (below + 0.5 * wn)).sum() / (wp.sum() * wn.sum()))


def load(pred_dir, task, train, cond, split, suffix, test):
    f = pred_dir / f"{task}__{train}__{cond}__split{split}{suffix}__on_{test}.pred.tsv"
    return pd.read_csv(f, sep="\t", usecols=["protein", "pos", "y", "y_pred"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pred_dir", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("-B", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tasks", nargs="*", default=TASKS)
    a = ap.parse_args()

    # the weighted estimator against sklearn, once, on real scores and random integer weights
    d = load(a.pred_dir, "methylation_r", "replica", "baseline", 0, "", "rebuilt")
    w = np.random.default_rng(1).integers(0, 4, len(d)).astype(float)
    ref = roc_auc_score(d.y, d.y_pred, sample_weight=w)
    assert abs(Scored(d.y.values, d.y_pred.values).auroc(w) - ref) < 1e-12, ref

    rows = []
    for task in a.tasks:
        for train, test in CELLS:
            t0 = time.time()
            base = {s: load(a.pred_dir, task, train, "baseline", s, "", test) for s in SPLITS}
            aug = {(f, s): load(a.pred_dir, task, train, "ppi", s, suf, test)
                   for f, suf in FAMILIES.items() for s in SPLITS}
            for (f, s), x in aug.items():
                b = base[s]
                assert (x.protein.values == b.protein.values).all() and (x.pos.values == b.pos.values).all() \
                    and (x.y.values == b.y.values).all(), (task, train, test, f, s)
            members = {s: set(base[s].protein) for s in SPLITS}
            prots = sorted(set().union(*members.values()))
            index = {p: i for i, p in enumerate(prots)}
            n = len(prots)
            shared = sum(1 for p in prots if sum(p in members[s] for s in SPLITS) > 1)
            pidx = {s: base[s].protein.map(index).to_numpy() for s in SPLITS}
            B0 = {s: Scored(base[s].y.values, base[s].y_pred.values) for s in SPLITS}
            A0 = {k: Scored(x.y.values, x.y_pred.values) for k, x in aug.items()}

            point = {}
            for f in FAMILIES:
                per = [roc_auc_score(aug[(f, s)].y, aug[(f, s)].y_pred) - roc_auc_score(base[s].y, base[s].y_pred)
                       for s in SPLITS]
                ones = {s: np.ones(len(base[s])) for s in SPLITS}
                chk = [A0[(f, s)].auroc(ones[s]) - B0[s].auroc(ones[s]) for s in SPLITS]
                assert np.allclose(per, chk, atol=1e-12), (task, train, test, f)
                point[f] = per

            tag = f"{a.seed}|{task}|{train}|{test}".encode()
            rng = np.random.default_rng(int.from_bytes(hashlib.sha256(tag).digest()[:8], "little"))
            stat = {f: np.empty(a.B) for f in FAMILIES}
            kept = 0
            for _ in range(a.B):
                mult = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
                ws = {s: mult[pidx[s]] for s in SPLITS}
                if any((ws[s] * B0[s].pos).sum() == 0 or (ws[s] * ~B0[s].pos).sum() == 0 for s in SPLITS):
                    continue  # a split left without one of the classes; not observed in practice
                b_auc = {s: B0[s].auroc(ws[s]) for s in SPLITS}
                for f in FAMILIES:
                    stat[f][kept] = np.mean([A0[(f, s)].auroc(ws[s]) - b_auc[s] for s in SPLITS])
                kept += 1
            for f in FAMILIES:
                lo, hi = np.percentile(stat[f][:kept], [2.5, 97.5])
                rows.append(dict(task=task, train=train, eval=test, family=f, d_auroc=np.mean(point[f]),
                                 auroc_lo=lo, auroc_hi=hi, d_split0=point[f][0], d_split1=point[f][1],
                                 d_split2=point[f][2], n_boot=kept, n_proteins=n, n_shared_proteins=shared,
                                 n_sites=sum(len(base[s]) for s in SPLITS)))
            print(f"{task:20s} {train}->{test}  {n} proteins ({shared} shared)  {kept} resamples  "
                  f"{time.time() - t0:5.0f} s", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep="\t", index=False, float_format="%.6f")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
