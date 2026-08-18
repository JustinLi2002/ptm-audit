#!/usr/bin/env python3
"""Supplementary Figure S2 -- false-negative sensitivity.

(a) Mean Delta AUROC against alpha, both feature families, methylation K and R
    collapsed into a single unit (n = 7) as in the main analysis.
(b) Per task, the interval between the value under complete-annotation labels
    (alpha = 1, the published cell) and under the maximally adversarial
    assumption (alpha = 0), ordered by the former.

Panel (b) replaces an earlier plan to plot imputed false-negative burden
against harm.  There is no monotone relationship between them (Spearman
+0.17 for the language model, -0.10 for the interaction embedding, n = 8),
so such a panel would suggest a trend the data do not support.

    python make_figure_s2.py --results ../results --out ../figures
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402

MERGE = {"methylation_k": "methyl", "methylation_r": "methyl"}
LABEL = {
    "phosphorylation_st": "Phospho S/T", "phosphorylation_y": "Phospho Y",
    "ubiquitination_k": "Ubiquitination K", "sumoylation_k": "Sumoylation K",
    "acetylation_k": "Acetylation K", "methylation_k": "Methylation K/R",
    "methylation_r": "Methylation R", "glycosylation_n": "N-glycosylation",
}
ARMS = {"esm": ("Language model", "#1f4e79"),
        "ppi": ("Interaction embedding", "#c55a11")}


def fold(df):
    """Collapse the two nested methylation tasks into one unit (n = 7)."""
    d = df.copy()
    d["unit"] = d["task"].map(MERGE).fillna(d["task"])
    return d


def load(results, arm):
    p = Path(results) / f"fn_sensitivity_{arm}.tsv"
    d = pd.read_csv(p, sep="\t")
    return d[d["alpha"].notna() & d.get("target_inflation",
                                        pd.Series(index=d.index)).isna()] \
        if "target_inflation" in d.columns else d


def main(results, out):
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.9),
                           gridspec_kw={"width_ratios": [1, 1.25]})

    # ---- (a) alpha curves -------------------------------------------------
    for arm, (name, colour) in ARMS.items():
        d = fold(load(results, arm))
        per_unit = d.groupby(["alpha", "unit"])["d_auroc"].mean().reset_index()
        curve = per_unit.groupby("alpha")["d_auroc"].mean()
        ax[0].plot(curve.index, curve.values, "o-", color=colour,
                   lw=1.8, ms=5, label=name)
        a1, a0 = curve.loc[1.0], curve.loc[0.0]
        ax[0].annotate(f"{a0 / a1:.0%} retained",
                       xy=(0.0, a0), xytext=(0.06, a0 + 0.012),
                       color=colour, fontsize=8)

    ax[0].axhline(0, color="0.35", lw=0.8)
    ax[0].set_xlabel(r"$\alpha$   (0 = all of the depth–rate relationship"
                     "\nis study bias;  1 = annotation complete)", fontsize=8)
    ax[0].set_ylabel(r"mean $\Delta$AUROC, natural negatives")
    ax[0].set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax[0].legend(frameon=False, fontsize=8, loc="lower left")
    ax[0].set_title("a", loc="left", fontweight="bold")

    # ---- (b) per-task interval -------------------------------------------
    d = load(results, "esm")
    piv = (d[d["alpha"].isin([0.0, 1.0])]
           .pivot_table(index="task", columns="alpha", values="d_auroc")
           .sort_values(1.0))
    ypos = np.arange(len(piv))
    for y, (task, row) in zip(ypos, piv.iterrows()):
        lo, hi = row[1.0], row[0.0]
        flips = lo < 0 <= hi
        ax[1].plot([lo, hi], [y, y], "-",
                   color="#b03a2e" if flips else "0.6", lw=2.4, zorder=1)
        ax[1].scatter([lo], [y], s=34, color=ARMS["esm"][1], zorder=3)
        ax[1].scatter([hi], [y], s=34, facecolor="white", zorder=3,
                      edgecolor="#b03a2e" if flips else ARMS["esm"][1])

    ax[1].axvline(0, color="0.35", lw=0.8)
    ax[1].set_yticks(ypos)
    ax[1].set_yticklabels([LABEL.get(t, t) for t in piv.index], fontsize=8)
    ax[1].set_xlabel(r"$\Delta$AUROC, language model channel")
    ax[1].set_title("b", loc="left", fontweight="bold")
    ax[1].scatter([], [], s=34, color=ARMS["esm"][1], label=r"$\alpha=1$")
    ax[1].scatter([], [], s=34, facecolor="white",
                  edgecolor=ARMS["esm"][1], label=r"$\alpha=0$")
    ax[1].legend(frameon=False, fontsize=8, loc="lower right")

    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"figure_s2_fn_sensitivity.{ext}", dpi=300)
    print(f"wrote {out}/figure_s2_fn_sensitivity.pdf")

    print("\nvalues used in panel a:")
    for arm in ARMS:
        d = fold(load(results, arm))
        c = (d.groupby(["alpha", "unit"])["d_auroc"].mean().reset_index()
             .groupby("alpha")["d_auroc"].mean())
        print(f"  {arm}: " + "  ".join(f"a={a:.2f} {v:+.4f}"
                                       for a, v in c.items()))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="../results")
    ap.add_argument("--out", default="../figures")
    a = ap.parse_args()
    main(a.results, a.out)
