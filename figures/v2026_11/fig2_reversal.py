"""Figure 2. The donor rule reverses the measured contribution of protein-level features."""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from figstyle import (CELLS, CH, CH_LABEL, GRID, INK, INK2, MM, MUTED, RES, SHORT, TASKS,
                      hgrid, label, save, setup, vgrid)

FAMS = ["interaction", "esm2", "prott5"]
CODE = {("replica", "replica"): "TT", ("replica", "rebuilt"): "TU",
        ("rebuilt", "replica"): "UT", ("rebuilt", "rebuilt"): "UU"}
TITLE = {"TT": "TT  threshold train · threshold test", "TU": "TU  threshold train · unrestricted test",
         "UT": "UT  unrestricted train · threshold test", "UU": "UU  unrestricted train · unrestricted test"}


def load():
    v = pd.read_csv(RES / "crosseval_verify.tsv", sep="\t")
    v["cell"] = [CODE[(a, b)] for a, b in zip(v.train, v["eval"])]
    return v


def worked_example(ax_roc, ax_prc, v):
    """Acetylation K, ESM-2: absolute values of baseline and augmented model per cell."""
    d = v[(v.task == "acetylation_k") & (v.family == "esm2")]
    g = d.groupby("cell")[["auroc_base", "auroc", "auprc_base", "auprc", "pos_rate"]].mean()
    order = ["TT", "TU", "UT", "UU"]
    y = np.arange(4)[::-1]
    for ax, (b, a, name) in [(ax_roc, ("auroc_base", "auroc", "AUROC")), (ax_prc, ("auprc_base", "auprc", "AUPRC"))]:
        for yy, c in zip(y, order):
            x0, x1 = g.loc[c, b], g.loc[c, a]
            ax.annotate("", xy=(x1, yy), xytext=(x0, yy),
                        arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.0, shrinkA=3, shrinkB=4,
                                        mutation_scale=7))
            ax.scatter(x0, yy, s=28, facecolor="white", edgecolor=INK2, linewidth=1.0, zorder=3)
            ax.scatter(x1, yy, s=28, color=CH["esm2"], edgecolor="white", linewidth=0.8, zorder=3)
        ax.set_yticks(y, order if ax is ax_roc else [""] * 4)
        ax.set_xlabel(name)
        ax.set_ylim(-0.6, 3.6)
        vgrid(ax)
    ax_roc.set_xlim(0.55, 0.92)
    ax_prc.set_xlim(0.15, 0.88)
    for yy, c in zip(y, order):
        ax_prc.text(1.02, yy, f"{100 * g.loc[c, 'pos_rate']:.0f}% positive", transform=ax_prc.get_yaxis_transform(),
                    fontsize=6, color=INK2, va="center")
    ax_roc.set_title("Acetylation K, ESM-2", loc="left", fontsize=7)
    h = [Line2D([], [], ls="", marker="o", mfc="white", mec=INK2, ms=4.5, label="sequence only"),
         Line2D([], [], ls="", marker="o", color=CH["esm2"], ms=4.5, label="+ ESM-2 vector")]
    ax_prc.legend(handles=h, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, handletextpad=0.3,
                  columnspacing=1.0, borderaxespad=0.2)


def cell_panel(ax, v, cell, show_x, ylim):
    d = v[v.cell == cell]
    x = np.arange(len(TASKS))
    off = {"interaction": -0.22, "esm2": 0.0, "prott5": 0.22}
    ax.axhline(0, color=INK2, lw=0.7, zorder=1)
    for f in FAMS:
        real = d[d.family == f].groupby("task").d_auroc
        perm = d[d.family == f + "_perm"].groupby("task").d_auroc.mean()
        for i, t in enumerate(TASKS):
            pts = d[(d.family == f) & (d.task == t)].d_auroc.values
            ax.scatter(np.full(len(pts), i + off[f]), pts, s=3, color=MUTED, zorder=2, linewidth=0)
            m = real.mean()[t]
            wins = m > perm[t]
            if wins:
                ax.scatter(i + off[f], m, s=24, color=CH[f], edgecolor="white", linewidth=0.7, zorder=3)
            else:
                ax.scatter(i + off[f], m, s=24, marker="^", facecolor="white", edgecolor=CH[f],
                           linewidth=1.1, zorder=3)
    mean = {f: d[d.family == f].d_auroc.mean() for f in FAMS}
    ax.set_title(TITLE[cell], loc="left", fontsize=7)
    ax.text(0.99, 0.04 if cell != "TU" else 0.94,
            "mean  " + "   ".join(f"{CH_LABEL[f]} {mean[f]:+.3f}" for f in FAMS),
            transform=ax.transAxes, ha="right", va="bottom" if cell != "TU" else "top", fontsize=5.6, color=INK2)
    ax.set_xticks(x, [SHORT[t] for t in TASKS] if show_x else [""] * len(TASKS), rotation=40, ha="right")
    ax.set_xlim(-0.6, len(TASKS) - 0.4)
    ax.set_ylim(*ylim)
    hgrid(ax)


def uu_auprc(ax, v):
    d = v[v.cell == "UU"]
    x = np.arange(len(TASKS))
    w = 0.26
    for k, f in enumerate(FAMS):
        m = d[d.family == f].groupby("task").d_auprc.mean().reindex(TASKS)
        ax.bar(x + (k - 1) * w, m.values, width=w * 0.88, color=CH[f], label=CH_LABEL[f], zorder=2)
    ax.axhline(0, color=INK2, lw=0.7)
    ax.set_xticks(x, [SHORT[t] for t in TASKS], rotation=40, ha="right")
    ax.set_ylabel("ΔAUPRC")
    ax.set_title("UU  unrestricted train · unrestricted test, ΔAUPRC", loc="left", fontsize=7)
    hgrid(ax)
    ax.legend(loc="upper right", ncol=3, handlelength=1.0, columnspacing=1.0)


def main():
    setup()
    v = load()
    fig = plt.figure(figsize=(180 * MM, 200 * MM))
    gs = fig.add_gridspec(4, 2, height_ratios=[0.8, 1, 1, 0.9], hspace=0.62, wspace=0.18)
    top = gs[0, :].subgridspec(1, 3, width_ratios=[1, 1, 0.55], wspace=0.12)
    a1 = fig.add_subplot(top[0, 0]); a2 = fig.add_subplot(top[0, 1])
    worked_example(a1, a2, v)
    label(a1, "a", x=-0.18, y=1.1)
    lo = min(v[v.family.isin(FAMS)].d_auroc.min(), -0.3) - 0.02
    hi = v[v.family.isin(FAMS)].d_auroc.max() + 0.03
    axes = {}
    for k, cell in enumerate(["TT", "TU", "UT", "UU"]):
        ax = fig.add_subplot(gs[1 + k // 2, k % 2])
        cell_panel(ax, v, cell, show_x=(k // 2 == 1), ylim=(lo, hi))
        if k % 2 == 1:
            ax.tick_params(labelleft=False)
        label(ax, "bcde"[k], x=-0.1 if k % 2 == 0 else -0.05, y=1.05)
        axes[cell] = ax
    fig.text(0.045, 0.53, "ΔAUROC from adding the protein vector", rotation=90, ha="center",
             va="center", fontsize=7)
    hl = [Line2D([], [], ls="", marker="o", color=CH[f], ms=4.5, label=CH_LABEL[f]) for f in FAMS]
    hl += [Line2D([], [], ls="", marker="o", color=INK2, ms=4.5, label="real beats permuted"),
           Line2D([], [], ls="", marker="^", mfc="white", mec=INK2, ms=4.5, label="permuted wins"),
           Line2D([], [], ls="", marker="o", color=MUTED, ms=2, label="single partition")]
    axes["TT"].legend(handles=hl, loc="lower left", bbox_to_anchor=(0, 1.13), ncol=6, handletextpad=0.2,
                      columnspacing=0.9, borderaxespad=0)
    f_ax = fig.add_subplot(gs[3, :])
    uu_auprc(f_ax, v)
    label(f_ax, "f", x=-0.065, y=1.05)
    save(fig, "Figure2_reversal")
    # numbers quoted in the caption
    g = v[v.family.isin(FAMS)].groupby(["family", "cell"])[["d_auroc", "d_auprc"]].mean().round(3)
    print(g)


if __name__ == "__main__":
    main()
