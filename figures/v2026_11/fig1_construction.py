"""Figure 1. Construction and coverage."""
import re

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

from figstyle import AXIS, INK, INK2, MUTED, MM, RES, SHORT, TASKS, vgrid, label, save, setup


def r0(v):
    """Percent to an integer string, halves rounded up (60.5 -> 61, as in the manuscript)."""
    return f"{int(v + 0.5)}"


# Table 1 of the manuscript: published pure-positive, threshold reconstruction, unrestricted
PUREPOS = {
    "phosphorylation_st": (4.4, 4.0, 0.0), "phosphorylation_y": (2.4, 3.6, 0.4),
    "acetylation_k": (22.5, 20.9, 0.1), "methylation_k": (25.2, 25.5, 0.1),
    "methylation_r": (29.8, 29.6, 0.2), "sumoylation_k": (17.9, 17.2, 0.0),
    "ubiquitination_k": (13.5, 12.5, 0.1), "glycosylation_n": (38.3, 37.5, 0.1)}


def coverage():
    """Share of candidate residues on proteins with merged depth < 10 (analysis/unobserved_weight.py)."""
    name = {"Phosphorylation S/T": "phosphorylation_st", "Phosphorylation Y": "phosphorylation_y",
            "Acetylation K": "acetylation_k", "Methylation K/R": "methylation_k",
            "Methylation R": "methylation_r", "Sumoylation K": "sumoylation_k",
            "Ubiquitination K": "ubiquitination_k", "N-Glycosylation N": "glycosylation_n"}
    out = {}
    for line in open(RES / "unobserved_weight.txt", encoding="utf-8"):
        m = re.match(r"(\S.*?)\s{2,}([\d.]+)\s+(\d+)\s+([\d.]+)\s+(\d+)\s*$", line)
        if m and m.group(1) in name:
            out[name[m.group(1)]] = (100 * float(m.group(4)), 100 * float(m.group(2)))
    assert len(out) == 8, out
    return out  # task: (whole reference proteome %, proteins with >=1 site of the task %)


def protein(ax, y, n_cand, pos, sampled_neg, title, note):
    """One protein drawn as a row of candidate residues."""
    xs = np.linspace(0.05, 0.95, n_cand)
    ax.plot([0.02, 0.98], [y, y], color=AXIS, lw=3, solid_capstyle="round", zorder=1)
    for i, x in enumerate(xs):
        if i in pos:
            ax.scatter(x, y, s=24, color=INK, zorder=3, edgecolor="white", linewidth=0.8)
        elif i in sampled_neg:
            ax.scatter(x, y, s=24, facecolor="white", edgecolor=INK, linewidth=0.9, zorder=3)
        else:
            ax.scatter(x, y, s=9, color=MUTED, marker="|", zorder=2, linewidth=0.8)
    ax.text(0.02, y + 0.06, title, fontsize=6.5, color=INK, va="bottom")
    ax.text(0.98, y + 0.06, note, fontsize=6.2, color=INK2, va="bottom", ha="right")


def schematic(ax):
    ax.set_xlim(0, 1); ax.set_ylim(-0.14, 1); ax.axis("off")
    n = 22
    pos_a = {1, 3, 4, 6, 9, 11, 13, 15, 18, 20}     # a donor protein: ten annotated sites
    pos_b = {3, 10, 16}                             # three annotated sites, below the threshold
    neg_a = {0, 5, 8, 12, 17, 21}
    neg_b = {0, 5, 7, 12, 14, 19, 21}
    ax.text(0.0, 0.97, "Threshold construction (≥10-site donor rule)", fontsize=7, fontweight="bold", va="top")
    protein(ax, 0.76, n, pos_a, neg_a, "Protein A, depth 10", "negatives sampled")
    protein(ax, 0.58, n, pos_b, set(), "Protein B, depth 3", "no negatives sampled")
    ax.text(0.5, 0.52, "every retained site of protein B is positive", fontsize=6.2, color=INK2,
            ha="center", va="top", style="italic")
    ax.text(0.0, 0.42, "Unrestricted construction", fontsize=7, fontweight="bold", va="top")
    protein(ax, 0.23, n, pos_a, neg_a, "Protein A", "negatives sampled")
    protein(ax, 0.05, n, pos_b, neg_b, "Protein B", "negatives sampled")
    ky = -0.10
    ax.scatter(0.03, ky, s=22, color=INK)
    ax.text(0.055, ky, "positive (annotated)", fontsize=6, va="center")
    ax.scatter(0.43, ky, s=22, facecolor="white", edgecolor=INK, linewidth=0.9)
    ax.text(0.455, ky, "sampled negative", fontsize=6, va="center")
    ax.scatter(0.76, ky, s=9, color=MUTED, marker="|")
    ax.text(0.78, ky, "not retained", fontsize=6, va="center")


def purepos(ax):
    y = np.arange(len(TASKS))[::-1]
    h = 0.26
    cols = [("Published", "#c3c2b7"), ("Threshold reconstruction", INK2),
            ("Unrestricted reconstruction", "#2a78d6")]
    for k, (lab, c) in enumerate(cols):
        v = [PUREPOS[t][k] for t in TASKS]
        ax.barh(y + (1 - k) * h, v, height=h * 0.86, color=c, label=lab, zorder=2)
    ax.text(41, y[0] + 0.1, "Unrestricted reconstruction:\n0.0–0.4% in every task", fontsize=6,
            color=INK2, ha="right", va="center")
    ax.set_yticks(y, [SHORT[t] for t in TASKS])
    ax.set_xlabel("Training sites on pure-positive proteins (%)")
    ax.set_xlim(0, 42)
    vgrid(ax)
    ax.legend(loc="upper left", bbox_to_anchor=(0, -0.2), handlelength=1.0, borderaxespad=0)


def coverage_panel(ax, cov):
    y = np.arange(len(TASKS))[::-1]
    h = 0.34
    for k, (lab, c) in enumerate([("All unannotated candidate sites in the reference proteome", INK2),
                                  ("Unannotated candidate sites on proteins with ≥1 site of the task", "#a8a7a0")]):
        v = [cov[t][k] for t in TASKS]
        ax.barh(y + (0.5 - k) * h, v, height=h * 0.86, color=c, label=lab, zorder=2)
        for yy, vv in zip(y, v):
            ax.text(vv + 1, yy + (0.5 - k) * h, r0(vv), fontsize=5.6, va="center", color=INK2)
    ax.set_yticks(y, [SHORT[t] for t in TASKS])
    ax.set_xlim(0, 112)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Unannotated candidate sites on proteins with merged depth < 10 (%)")
    vgrid(ax)
    ax.legend(loc="upper left", bbox_to_anchor=(0, -0.2), handlelength=1.0, borderaxespad=0)


def design(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    x0, y0, w, h = 0.36, 0.20, 0.31, 0.28
    names = {(0, 0): "TT", (0, 1): "TU", (1, 0): "UT", (1, 1): "UU"}
    for r in range(2):
        for c in range(2):
            ax.add_patch(Rectangle((x0 + c * w, y0 + (1 - r) * h), w, h,
                                   facecolor="#e7e6e1" if (r, c) == (0, 0) else "#f4f3f0",
                                   edgecolor="white", linewidth=2))
            ax.text(x0 + c * w + w / 2, y0 + (1 - r) * h + h / 2, names[(r, c)], ha="center",
                    va="center", fontsize=9, fontweight="bold", color=INK)
    ax.text(x0 + w / 2, y0 + 2 * h + 0.03, "threshold", ha="center", fontsize=6.5)
    ax.text(x0 + 1.5 * w, y0 + 2 * h + 0.03, "unrestricted", ha="center", fontsize=6.5)
    ax.text(x0 + w, y0 + 2 * h + 0.11, "Evaluation construction", ha="center", fontsize=7, fontweight="bold")
    ax.text(x0 - 0.02, y0 + 1.5 * h, "threshold", ha="right", va="center", fontsize=6.5)
    ax.text(x0 - 0.02, y0 + 0.5 * h, "unrestricted", ha="right", va="center", fontsize=6.5)
    ax.text(0.03, y0 + h, "Training construction", ha="center", va="center", fontsize=7,
            fontweight="bold", rotation=90)
    ax.text(x0 + w, y0 - 0.05, "Same three protein-disjoint partitions in every cell;\n"
                               "each trained ensemble is scored on both test constructions.",
            fontsize=6.2, color=INK2, va="top", ha="center")


def main():
    setup()
    cov = coverage()
    fig = plt.figure(figsize=(180 * MM, 132 * MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], height_ratios=[1, 1.05], hspace=0.55, wspace=0.42)
    a = fig.add_subplot(gs[0, 0]); schematic(a); label(a, "a", x=-0.02, y=1.0)
    b = fig.add_subplot(gs[0, 1]); purepos(b); label(b, "b", x=-0.30)
    c = fig.add_subplot(gs[1, 0]); coverage_panel(c, cov); label(c, "c", x=-0.22)
    d = fig.add_subplot(gs[1, 1]); design(d); label(d, "d", x=-0.02, y=1.0)
    save(fig, "Figure1_construction_coverage")
    print({t: tuple(r0(v) for v in cov[t]) for t in TASKS})


if __name__ == "__main__":
    main()
