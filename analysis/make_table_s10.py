#!/usr/bin/env python3
"""Supplementary Table S10 -- false-negative sensitivity, per task.

Reports Delta AUROC under the published labels (alpha = 1), under two uniform
assumptions about the fraction of true sites that went unannotated, and at the
ceiling of what the parameterisation can express for each task (alpha = 0).

A target inflation is only comparable across tasks where every task can reach
it.  The ceiling is set by the dynamic range of the observed positive rate and
differs by task, so targets above the smallest ceiling are marked.

    python make_table_s10.py --results ../results --out ../results
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

LABEL = {
    "phosphorylation_st": "Phosphorylation S/T", "phosphorylation_y": "Phosphorylation Y",
    "ubiquitination_k": "Ubiquitination K", "sumoylation_k": "Sumoylation K",
    "acetylation_k": "Acetylation K", "methylation_k": "Methylation K",
    "methylation_r": "Methylation R", "glycosylation_n": "N-glycosylation",
}
FOLD = {"methylation_k": "methyl", "methylation_r": "methyl"}


def build(results, arm):
    d = pd.read_csv(Path(results) / f"fn_sensitivity_{arm}.tsv", sep="\t")
    grid = d[d["mode"] == "alpha"]
    infl = d[d["mode"] == "inflation"]

    out = pd.DataFrame(index=sorted(d["task"].unique()))
    out["ceiling"] = d.groupby("task")["max_inflation"].mean()
    out["burden_a0"] = (grid[grid.alpha == 0].groupby("task")
                        .apply(lambda g: (g.exp_flipped / g.n_sites).mean(),
                               include_groups=False))
    out["a1"] = grid[grid.alpha == 1].groupby("task")["d_auroc"].mean()

    for lam in sorted(infl["target_inflation"].unique()):
        sub = infl[infl.target_inflation == lam]
        out[f"lam{lam:g}"] = sub.groupby("task")["d_auroc"].mean()
        reached = sub.groupby("task")["inflation_reached"].all()
        out[f"lam{lam:g}_capped"] = ~reached
        out[f"lam{lam:g}_actual"] = sub.groupby("task")["actual_inflation"].mean()

    out["a0"] = grid[grid.alpha == 0].groupby("task")["d_auroc"].mean()
    # A ratio is only meaningful where the denominator is an effect and the
    # effect keeps its sign.  Below this magnitude the "retained fraction" is
    # dominated by its denominator and reports numbers like -3900%; where the
    # effect crosses zero a negative "retained" percentage is meaningless.
    out["retained_a0"] = ((out["a0"] / out["a1"])
                          .where((out["a1"] < -0.01) & (out["a0"] < 0)))
    return out.sort_values("a1")


def build_both(results):
    """Common row order across arms so the two tables can be read together."""
    tables = {arm: build(results, arm) for arm in ("esm", "ppi")}
    order = list(tables["esm"].index)
    return {arm: t.loc[order] for arm, t in tables.items()}


def render(t, arm):
    lam_cols = [c for c in t.columns
                if c.startswith("lam") and not c.endswith(("_capped", "_actual"))]
    lines = [f"### {arm}", "",
             "| PTM task | ceiling | burden at α=0 | ΔAUROC α=1 | "
             + " | ".join(f"λ={c[3:]}" for c in lam_cols)
             + " | ΔAUROC α=0 | retained |",
             "|---|---|---|---|" + "---|" * (len(lam_cols) + 2)]
    for task, r in t.iterrows():
        cells = []
        for c in lam_cols:
            mark = "\u2020" if r[f"{c}_capped"] else ""
            cells.append(f"{r[c]:+.4f}{mark}")
        ret = "—" if pd.isna(r["retained_a0"]) else f"{r['retained_a0']:.0%}"
        lines.append(
            f"| {LABEL.get(task, task)} | {r['ceiling']:.2f}× | "
            f"{r['burden_a0']:.1%} | {r['a1']:+.4f} | "
            + " | ".join(cells) + f" | {r['a0']:+.4f} | {ret} |")

    # 折叠后的合计,与正文一致
    t2 = t.copy()
    t2["unit"] = [FOLD.get(i, i) for i in t2.index]
    u = t2.groupby("unit")[["a1", "a0"] + lam_cols].mean()
    lines += ["", f"Mean over 7 units (methylation K and R collapsed): "
                  f"α=1 {u['a1'].mean():+.4f}, "
              + ", ".join(f"λ={c[3:]} {u[c].mean():+.4f} "
                          f"({u[c].mean() / u['a1'].mean():.0%})"
                          for c in lam_cols)
              + f", α=0 {u['a0'].mean():+.4f} "
                f"({u['a0'].mean() / u['a1'].mean():.0%}).",
              "",
              "† target above this task's ceiling; the value shown is at α = 0, "
              "where the achieved inflation is the ceiling rather than the "
              "target.", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="../results")
    ap.add_argument("--out", default="../results")
    a = ap.parse_args()

    md = ["# Supplementary Table S10", "",
          "Delta AUROC of the added channel under natural negatives, as a "
          "function of assumed annotation incompleteness. Threshold-sampled "
          "training, protein-disjoint evaluation, mean over three partitions.",
          ""]
    tables = build_both(a.results)
    for arm, name in (("esm", "Language model channel"),
                      ("ppi", "Interaction embedding channel")):
        t = tables[arm]
        t.to_csv(Path(a.out) / f"table_s10_{arm}.tsv", sep="\t")
        md.append(render(t, name))
        print(render(t, name))

    p = Path(a.out) / "table_s10.md"
    p.write_text("\n".join(md))
    print(f"wrote {p}")
