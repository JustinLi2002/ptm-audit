"""Figure 2 (v7). As fig2_reversal_v6.py, with larger in-panel numbers and v7 terminology.

g, h: the mean of each training construction is written next to its points at about 7 pt final size, in place of
the row of means at the panel foot. b-e: the row of means is enlarged and uses a minus sign. 'single partition'
becomes 'single split' (a split is one of the three protein-disjoint replicates)."""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from fig2_reversal import FAMS, cell_panel, load, uu_auprc, worked_example
from figstyle import CH, CH_LABEL, INK, INK2, MM, MUTED, hgrid, label, save, setup

RUNS = Path(__file__).resolve().parents[2] / "ratio_matched" / "real" / "runs_real.tsv"
# construction, tick label, marker, face colour (edge is always the ESM-2 colour)
GROUPS = [("T", "T", "o", CH["esm2"]), ("matched", "Ratio-matched", "s", CH["esm2"]), ("U", "U", "o", "white")]
# the figure is inserted at 15.0 cm of a 15.7 cm drawing, so 7.3 pt here prints at about 7 pt
NUM_PT = 7.3
MEAN_ROW_PT = 6.8


def minus(s):
    return s.replace("-", "−")


def matched_panel(ax, r, task, title):
    d = r[r.task == task]
    p = d.pivot(index="split", columns="construction", values="dauroc_unrestricted")
    assert p.shape == (3, 3) and not p.isna().any().any()
    rate = d.groupby("construction").pos_train.mean()
    x = np.arange(len(GROUPS))
    ax.axhline(0, color=INK2, lw=0.7, zorder=1)
    for s in p.index:
        ax.plot(x, [p.loc[s, g] for g, *_ in GROUPS], color=MUTED, lw=0.6, alpha=0.8, zorder=2)
    lo, hi = min(p.values.min(), 0), max(p.values.max(), 0)
    pad = 0.1 * (hi - lo)
    for xi, (g, _, m, fc) in zip(x, GROUPS):
        ax.scatter(np.full(len(p), xi), p[g], s=20, marker=m, facecolor=fc, edgecolor=CH["esm2"],
                   linewidth=1.0, zorder=3)
        mean = p[g].mean()
        ax.plot([xi - 0.18, xi + 0.18], [mean] * 2, color=INK, lw=1.4, zorder=4)
        # the mean goes on the far side of the points from zero, away from the lines joining the splits
        below = mean < 0
        ax.text(xi, (p[g].min() - 0.4 * pad) if below else (p[g].max() + 0.4 * pad), minus(f"{mean:+.3f}"),
                ha="center", va="top" if below else "bottom", fontsize=NUM_PT, color=INK, zorder=5)
    ax.set_xticks(x, [f"{lab}\n{100 * rate[g]:.0f}% positive" for g, lab, *_ in GROUPS])
    ax.set_xlim(-0.5, len(GROUPS) - 0.5)
    ax.set_ylim(lo - 3.0 * pad, hi + 3.0 * pad)
    ax.set_xlabel("Training construction (training positive rate)")
    ax.set_ylabel("ΔAUROC from ESM-2\n(unrestricted test)")
    ax.set_title(title, loc="left", fontsize=7)
    hgrid(ax)
    print(task, "\n", p.round(4).to_string(), "\n means", p.mean().round(4).to_dict())


def main():
    setup()
    v = load()
    r = pd.read_csv(RUNS, sep="\t")
    fig = plt.figure(figsize=(180 * MM, 222 * MM))
    # panel rows (a, b-c, d-e, f, g-h) alternate with spacer rows sized in mm, so that rotated task labels and
    # panel titles keep a fixed clearance however many rows the figure has
    ROWS = [22, 17, 28, 12, 28, 21, 25, 21, 30]
    gs = fig.add_gridspec(len(ROWS), 2, height_ratios=ROWS, hspace=0, wspace=0.18, top=0.97, bottom=0.03)
    top = gs[0, :].subgridspec(1, 3, width_ratios=[1, 1, 0.55], wspace=0.12)
    a1 = fig.add_subplot(top[0, 0]); a2 = fig.add_subplot(top[0, 1])
    worked_example(a1, a2, v)
    label(a1, "a", x=-0.18, y=1.1)
    lo = min(v[v.family.isin(FAMS)].d_auroc.min(), -0.3) - 0.02
    hi = v[v.family.isin(FAMS)].d_auroc.max() + 0.03
    axes = {}
    for k, cell in enumerate(["TT", "TU", "UT", "UU"]):
        ax = fig.add_subplot(gs[2 + 2 * (k // 2), k % 2])
        cell_panel(ax, v, cell, show_x=(k // 2 == 1), ylim=(lo, hi))
        for t in ax.texts:
            if t.get_text().startswith("mean"):
                t.set_text(minus(t.get_text()))
                t.set_fontsize(MEAN_ROW_PT)
        if k % 2 == 1:
            ax.tick_params(labelleft=False)
        label(ax, "bcde"[k], x=-0.1 if k % 2 == 0 else -0.05, y=1.05)
        axes[cell] = ax
    yb = (axes["TT"].get_position().y1 + axes["UT"].get_position().y0) / 2
    fig.text(0.045, yb, "ΔAUROC from adding the protein vector", rotation=90, ha="center", va="center",
             fontsize=7)
    hl = [Line2D([], [], ls="", marker="o", color=CH[f], ms=4.5, label=CH_LABEL[f]) for f in FAMS]
    hl += [Line2D([], [], ls="", marker="o", color=INK2, ms=4.5, label="real beats permuted"),
           Line2D([], [], ls="", marker="^", mfc="white", mec=INK2, ms=4.5, label="permuted wins"),
           Line2D([], [], ls="", marker="o", color=MUTED, ms=2, label="single split")]
    axes["TT"].legend(handles=hl, loc="lower left", bbox_to_anchor=(0, 1.13), ncol=6, handletextpad=0.2,
                      columnspacing=0.9, borderaxespad=0)
    f_ax = fig.add_subplot(gs[6, :])
    uu_auprc(f_ax, v)
    label(f_ax, "f", x=-0.065, y=1.05)
    gh = gs[8, :].subgridspec(1, 2, wspace=0.36)
    for k, (task, title) in enumerate([("acetylation_k", "Acetylation K, ESM-2: ratio-matched control"),
                                       ("phosphorylation_st", "Phosphorylation S/T, ESM-2: ratio-matched control")]):
        ax = fig.add_subplot(gh[0, k])
        matched_panel(ax, r, task, title)
        label(ax, "gh"[k], x=-0.24, y=1.05)
    save(fig, "Figure2_reversal_v7")


if __name__ == "__main__":
    main()
