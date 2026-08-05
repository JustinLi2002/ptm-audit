#!/usr/bin/env python3
"""
make_figures.py — Figures 2, 3 and 5.

All values are hard-coded from the manuscript tables, so this runs anywhere.
Figure 4 needs per-site predictions and is generated separately on the cluster.

Outputs fig2.png/.pdf, fig3.png/.pdf, fig5.png/.pdf at 300 dpi.
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.8, "axes.edgecolor": "#444444",
    "xtick.color": "#444444", "ytick.color": "#444444",
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "axes.labelsize": 10, "legend.fontsize": 8.5,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
})

OUTDIR = os.path.dirname(os.path.abspath(__file__)) + "/figures"
os.makedirs(OUTDIR, exist_ok=True)

BLUE, GREY, ORANGE, PINK = "#1F6FB2", "#8A8A8A", "#D2691E", "#B5658C"
TASKS = ["Phospho S/T", "Phospho Y", "Acetylation K", "Methylation K/R",
         "Methylation R", "Sumoylation K", "Ubiquitination K", "N-Glyc N"]


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUTDIR}/{name}.{ext}")
    plt.close(fig)
    print(f"wrote {name}.png / .pdf")


# ─────────────────────────── Figure 2 ───────────────────────────
def figure2():
    id_full = [0.7729, 0.8020, 0.8765, 0.9571, 0.9450, 0.8201, 0.7990, 0.9247]
    id_rest = [0.7300, 0.7527, 0.7494, 0.7830, 0.8001, 0.7092, 0.7229, 0.5914]
    sq_full = [0.9510, 0.8710, 0.9049, 0.9503, 0.9306, 0.8615, 0.8814, 0.9868]
    sq_rest = [0.9524, 0.8696, 0.8847, 0.8854, 0.9233, 0.8370, 0.8712, 0.9832]

    fig, ax = plt.subplots(figsize=(7.2, 4.1))
    x = np.arange(len(TASKS))
    dx = 0.17

    for i in x:
        ax.plot([i - dx] * 2, [id_full[i], id_rest[i]], color=BLUE, lw=1.6, zorder=1)
        ax.plot([i + dx] * 2, [sq_full[i], sq_rest[i]], color=GREY, lw=1.6, zorder=1)
        ax.axvline(i + 0.5, color="#EEEEEE", lw=0.7, zorder=0)

    ax.scatter(x - dx, id_full, s=46, color=BLUE, zorder=3, label="identity, full test set")
    ax.scatter(x - dx, id_rest, s=42, facecolors="white", edgecolors=BLUE,
               linewidths=1.5, zorder=3, label="identity, restricted")
    ax.scatter(x + dx, sq_full, s=46, color=GREY, zorder=3, label="sequence model, full test set")
    ax.scatter(x + dx, sq_rest, s=42, facecolors="white", edgecolors=GREY,
               linewidths=1.5, zorder=3, label="sequence model, restricted")

    # annotate the N-glycosylation contrast, placed inside the axes
    ax.annotate("\u22120.333", xy=(7 - dx, (id_full[7] + id_rest[7]) / 2),
                xytext=(6.15, 0.665), color=BLUE, fontsize=9, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.1))
    ax.annotate("+0.0036", xy=(7 + dx, sq_rest[7]), xytext=(5.95, 0.945),
                color="#555555", fontsize=9, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="#777777", lw=1.1))

    ax.set_xticks(x)
    ax.set_xticklabels(TASKS, rotation=28, ha="right")
    ax.set_ylabel("AUROC")
    ax.set_ylim(0.55, 1.02)
    ax.set_xlim(-0.6, len(TASKS) - 0.4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower left", frameon=True, framealpha=0.95, edgecolor="#CCCCCC",
              ncol=2, columnspacing=1.1, handletextpad=0.5)
    save(fig, "fig2")


# ─────────────────────────── Figure 3 ───────────────────────────
def figure3():
    tasks = ["Acetylation K", "Methylation K/R", "Methylation R",
             "Sumoylation K", "Ubiquitination K", "N-Glyc N"]
    base = [0.9049, 0.9503, 0.9306, 0.8615, 0.8814, 0.9868]
    real = [0.9615, 0.9814, 0.9810, 0.9256, 0.9296, 0.9968]
    perm = [0.9606, 0.9821, 0.9813, 0.9242, 0.9277, 0.9970]
    cols = ["#1F6FB2", "#D2691E", "#B5658C", "#5B9BD5", "#7A4B2A", "#3D3D8F"]

    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    for i, t in enumerate(tasks):
        y = [base[i], real[i], perm[i]]
        ax.plot([0, 1, 2], y, "-o", color=cols[i], lw=1.5, ms=5, zorder=3)
        ax.text(2.06, perm[i], "  " + t, color=cols[i], fontsize=8.5,
                va="center", ha="left")

    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["sequence\nonly", "kinase\nfeature", "permuted\nkinase"])
    ax.set_xlim(-0.25, 3.05)
    ax.set_ylim(0.828, 1.005)
    ax.set_ylabel("AUROC")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#EEEEEE", lw=0.7)
    ax.set_axisbelow(True)
    save(fig, "fig3")


# ─────────────────────────── Figure 5 ───────────────────────────
def figure5():
    earlier = [0.0144, 0.0510, 0.0511, 0.0282, 0.0432, 0.0595, 0.0433, 0.0086]
    current = [-0.0028, 0.0040, -0.1067, -0.0596, -0.0453, -0.0304, -0.0079, -0.0029]
    rebuilt = [0.0136, 0.0136, 0.0471, 0.0112, 0.0181, 0.0299, 0.0274, 0.0038]
    parts = [[0.0136, 0.0148, 0.0124], [0.0161, 0.0131, 0.0119],
             [0.0436, 0.0519, 0.0458], [0.0127, 0.0093, 0.0116],
             [0.0110, 0.0254, 0.0180], [0.0233, 0.0282, 0.0380],
             [0.0290, 0.0276, 0.0252], [0.0051, 0.0039, 0.0025]]

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    y = np.arange(len(TASKS))[::-1]          # first task at the top
    off = 0.24

    ax.axvspan(-0.125, 0, color="#F6F1EE", zorder=0)
    ax.axvline(0, color="#222222", lw=1.3, zorder=2)
    for i in y:
        ax.axhline(i, color="#F0F0F0", lw=0.7, zorder=0)

    for k in range(len(TASKS)):
        yy = y[k]
        ax.scatter(parts[k], [yy - off] * 3, s=16, color=PINK, alpha=0.45, zorder=3)

    ax.scatter(earlier, y + off, s=52, color=BLUE, marker="o", zorder=4,
               label="earlier release, site-level split")
    ax.scatter(current, y, s=58, color=ORANGE, marker="^", zorder=4,
               label="current release, protein-disjoint split")
    ax.scatter(rebuilt, y - off, s=54, facecolors="white", edgecolors=PINK,
               linewidths=1.6, marker="s", zorder=5,
               label="threshold removed, protein-disjoint split")

    ax.set_yticks(y)
    ax.set_yticklabels(TASKS)
    ax.set_ylim(-0.75, len(TASKS) - 0.25)
    ax.set_xlim(-0.125, 0.078)
    ax.set_xlabel("\u0394 AUROC from adding the interaction channel")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="upper left", frameon=True, framealpha=0.96,
              edgecolor="#CCCCCC", handletextpad=0.6, borderpad=0.6)
    save(fig, "fig5")


if __name__ == "__main__":
    figure2()
    figure3()
    # figure5() superseded by analysis/make_figures_v2.py
