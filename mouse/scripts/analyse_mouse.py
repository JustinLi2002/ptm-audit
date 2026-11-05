#!/usr/bin/env python3
"""Mouse replication: outcomes and decision rules of PREREG_mouse.md, written before any result was seen.

Primary statistic: ESM-2 ΔAUROC (real minus sequence-only baseline), using the ensemble of the two
initialisations, averaged over the three splits.

The 95% interval comes from the joint protein bootstrap of bootstrap_20261102/crosseval_bootstrap_v2.py:
- 2,000 resamples;
- the same per-cell seeding;
- proteins drawn from the union of the three test partitions.

Rows A-D are read from the TT, TU and UU intervals. The matched control is read with the rule of the real-data
control (PREREG_real.md).

Usage: python analyse_mouse.py FETCH_DIR OUT_DIR
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

FETCH, OUT = Path(sys.argv[1]), Path(sys.argv[2])
RUNS = FETCH / "pdisjoint_runs_v2"
spec = importlib.util.spec_from_file_location(
    "bv2", Path(__file__).resolve().parents[2] / "analysis" / "crosseval_bootstrap_v2.py")
bv2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bv2)

TASK, SPLITS, B, SEED = "phosphorylation_st", (0, 1, 2), 2000, 0
NAME = {("replica", "replica"): "TT", ("replica", "rebuilt"): "TU", ("rebuilt", "replica"): "UT",
        ("rebuilt", "rebuilt"): "UU", ("matched", "replica"): "MT", ("matched", "rebuilt"): "MU"}
FAM = {"esm2": "ppi", "esm2_perm": "shuffled"}


def pred(train, cond, s, test):
    suf = "" if cond == "baseline" else "__esm"
    return pd.read_csv(RUNS / f"{TASK}__{train}__{cond}__split{s}{suf}__on_{test}.pred.tsv", sep="\t",
                       usecols=["protein", "pos", "y", "y_pred"])


rows = []
for (train, test), cell in NAME.items():
    base = {s: pred(train, "baseline", s, test) for s in SPLITS}
    aug = {(f, s): pred(train, c, s, test) for f, c in FAM.items() for s in SPLITS}
    for (f, s), x in aug.items():
        b = base[s]
        assert (x.protein.values == b.protein.values).all() and (x.pos.values == b.pos.values).all() \
            and (x.y.values == b.y.values).all(), (cell, f, s)
    members = {s: set(base[s].protein) for s in SPLITS}
    prots = sorted(set().union(*members.values()))
    index = {p: i for i, p in enumerate(prots)}
    n = len(prots)
    pidx = {s: base[s].protein.map(index).to_numpy() for s in SPLITS}
    B0 = {s: bv2.Scored(base[s].y.values, base[s].y_pred.values) for s in SPLITS}
    A0 = {k: bv2.Scored(x.y.values, x.y_pred.values) for k, x in aug.items()}
    point = {f: [roc_auc_score(aug[(f, s)].y, aug[(f, s)].y_pred) - roc_auc_score(base[s].y, base[s].y_pred)
                 for s in SPLITS] for f in FAM}
    prc = {f: [average_precision_score(aug[(f, s)].y, aug[(f, s)].y_pred)
               - average_precision_score(base[s].y, base[s].y_pred) for s in SPLITS] for f in FAM}
    real_gt_perm = sum(roc_auc_score(aug[("esm2", s)].y, aug[("esm2", s)].y_pred)
                       > roc_auc_score(aug[("esm2_perm", s)].y, aug[("esm2_perm", s)].y_pred) for s in SPLITS)

    tag = f"{SEED}|{TASK}|{train}|{test}".encode()
    rng = np.random.default_rng(int.from_bytes(hashlib.sha256(tag).digest()[:8], "little"))
    stat = {f: np.empty(B) for f in FAM}
    kept = 0
    for _ in range(B):
        mult = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
        ws = {s: mult[pidx[s]] for s in SPLITS}
        if any((ws[s] * B0[s].pos).sum() == 0 or (ws[s] * ~B0[s].pos).sum() == 0 for s in SPLITS):
            continue
        b_auc = {s: B0[s].auroc(ws[s]) for s in SPLITS}
        for f in FAM:
            stat[f][kept] = np.mean([A0[(f, s)].auroc(ws[s]) - b_auc[s] for s in SPLITS])
        kept += 1
    for f in FAM:
        lo, hi = np.percentile(stat[f][:kept], [2.5, 97.5])
        rows.append(dict(cell=cell, train=train, test=test, family=f, d_auroc=np.mean(point[f]), lo=lo, hi=hi,
                         d_split0=point[f][0], d_split1=point[f][1], d_split2=point[f][2],
                         splits_positive=sum(v > 0 for v in point[f]), d_auprc=np.mean(prc[f]),
                         real_gt_perm=real_gt_perm if f == "esm2" else np.nan, n_boot=kept, n_proteins=n,
                         base_auroc=np.mean([roc_auc_score(base[s].y, base[s].y_pred) for s in SPLITS]),
                         test_pos_rate=np.mean([base[s].y.mean() for s in SPLITS])))
    print(f"{cell}: done ({kept} resamples, {n} proteins)", flush=True)

r = pd.DataFrame(rows)
OUT.mkdir(parents=True, exist_ok=True)
r.to_csv(OUT / "mouse_results.tsv", sep="\t", index=False, float_format="%.6f")
e = r[r.family == "esm2"].set_index("cell")
above = lambda c: e.loc[c, "lo"] > 0
below = lambda c: e.loc[c, "hi"] < 0

if above("TT") and below("TU"):
    row = "A: direct replication: the feature helps on the threshold-sampled test and harms on the unrestricted test"
elif below("TU") and above("UU"):
    row = "B: training-construction effect replicated; the threshold-sampled test does not hide it"
elif below("TU"):
    row = "C: harm on the unrestricted test after threshold training, without the reversal"
else:
    row = "D: no harm; not replicated"

# matched, read with the rule of PREREG_real.md on the U test (MU against TU and UU), per split
mu = e.loc["MU", ["d_split0", "d_split1", "d_split2"]].to_numpy(float)
mT, mU, mM = e.loc["TU", "d_auroc"], e.loc["UU", "d_auroc"], e.loc["MU", "d_auroc"]
if not (mT < 0 and mU > 0):
    mrow = "1: T not negative or U not positive on the U test; matched reported descriptively"
elif (mu > 0).all() and abs(mM - mU) < abs(mM - mT):
    mrow = "2: class ratio and training-set size alone do not explain the contrast"
elif mM <= 0 or abs(mM - mT) < abs(mM - mU):
    mrow = "3: the mechanism attribution must be narrowed"
else:
    mrow = "4: mixed; effect sizes reported"

f = lambda x: f"{x:+.4f}".replace("-", "−")
lines = ["Mouse phosphorylation S/T, ESM-2 (PREREG_mouse.md). ΔAUROC = real − baseline, ensemble of two "
         "initialisations, mean of three splits; 95% joint protein bootstrap, 2,000 resamples.", ""]
for cell in NAME.values():
    x, p = e.loc[cell], r[(r.cell == cell) & (r.family == "esm2_perm")].iloc[0]
    lines.append(f"{cell}: ΔAUROC {f(x.d_auroc)} ({f(x.lo)}, {f(x.hi)}); splits > 0: {x.splits_positive}/3 "
                 f"[{', '.join(f(v) for v in (x.d_split0, x.d_split1, x.d_split2))}]; ΔAUPRC {f(x.d_auprc)}; "
                 f"permuted {f(p.d_auroc)}; real > permuted {int(x.real_gt_perm)}/3; baseline AUROC "
                 f"{x.base_auroc:.3f}; test positive rate {x.test_pos_rate:.3f}")
lines += ["", "Decision row: " + row, "Matched (secondary): row " + mrow]
ms = FETCH / "matched" / "matched_summary.tsv"
if ms.exists():
    m = pd.read_csv(ms, sep="\t").query("part == 'train'")
    lines.append("Pure-positive share of training sites, mean over splits: "
                 + ", ".join(f"{k} {100 * m[f'purepos_{k}'].mean():.1f}%" for k in ("T", "matched", "U"))
                 + f"; training positive rate T {m.rate_T.mean():.3f}, matched {m.rate_matched.mean():.3f}, "
                   f"U {m.rate_U.mean():.3f}")
(OUT / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("\n".join(lines))
