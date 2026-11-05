#!/usr/bin/env python3
"""Table 1 of the manuscript, from the deposited tables.

Panel a: results/crosseval_verify.tsv and results/crosseval_ci_v2.tsv (eight human tasks, three protein-level
features, four cells). Panel b: mouse/results/mouse_results.tsv (mouse phosphorylation S/T, ESM-2, the four cells read
by the rules fixed in mouse/PREREG_mouse.md).

    python analysis/table1_cells.py          # prints both panels and writes results/table1_cells.tsv
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FEATS = [("interaction", "Interaction"), ("esm2", "ESM-2"), ("prott5", "ProtT5")]
CELLS = [("replica", "replica", "TT"), ("replica", "rebuilt", "TU"), ("rebuilt", "replica", "UT"),
         ("rebuilt", "rebuilt", "UU")]


def f3(x):
    return f"{x:+.3f}".replace("-", "−")


v = pd.read_csv(ROOT / "results" / "crosseval_verify.tsv", sep="\t")
ci = pd.read_csv(ROOT / "results" / "crosseval_ci_v2.tsv", sep="\t")
rows = []
for fam, name in FEATS:
    for tr, ev, cell in CELLS:
        x = v[(v.family == fam) & (v.train == tr) & (v["eval"] == ev)].set_index(["task", "split"])
        p = v[(v.family == fam + "_perm") & (v.train == tr) & (v["eval"] == ev)].set_index(["task", "split"])
        c = ci[(ci.family == fam) & (ci.train == tr) & (ci["eval"] == ev)]
        tm = x.groupby("task")[["d_auroc", "d_auprc"]].mean()
        assert len(tm) == 8 and len(x) == 24 and len(c) == 8
        rows.append(dict(panel="a", feature=name, cell=cell, d_auroc=f3(tm.d_auroc.mean()),
                         d_auprc=f3(tm.d_auprc.mean()), tasks_with_gain_of_8=int((tm.d_auroc > 0).sum()),
                         task_intervals_above=int((c.auroc_lo > 0).sum()),
                         task_intervals_below=int((c.auroc_hi < 0).sum()), pairs_with_gain_of_24=int((x.d_auroc > 0).sum()),
                         real_gt_permuted_of_24=int(((x.auroc - p.auroc) > 0).sum())))
m = pd.read_csv(ROOT / "mouse" / "results" / "mouse_results.tsv", sep="\t")
m = m[m.family == "esm2"].set_index("cell")
for cell in ("TT", "TU", "UU", "MU"):
    r = m.loc[cell]
    rows.append(dict(panel="b", feature="ESM-2", cell=cell, d_auroc=f3(r.d_auroc), d_auprc=f3(r.d_auprc),
                     interval_95=f"{f3(r.lo)}, {f3(r.hi)}", splits_with_gain_of_3=int(r.splits_positive),
                     real_gt_permuted_of_3=int(r.real_gt_perm)))
t = pd.DataFrame(rows).convert_dtypes()
t.to_csv(ROOT / "results" / "table1_cells.tsv", sep="\t", index=False, lineterminator="\n")
pd.set_option("display.width", 220, "display.max_columns", 20)
print(t.astype("string").fillna("").to_string(index=False))
