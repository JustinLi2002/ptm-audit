"""Figure 3. Observations consistent with a sampling-induced depth association."""
from itertools import permutations

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from figstyle import (AXIS, CH, CH_LABEL, INK, INK2, MM, MUTED, RES, S3, SHORT, TASKS,
                      hgrid, label, save, setup, vgrid)

FAM = {"ppi": "interaction", "esm": "esm2"}
GROUP = [("phospho", "Phospho S/T/Y"), ("acet", "Acetyl K"), ("methyl", "Methyl K/R"),
         ("sumo", "Sumo K"), ("ubi", "Ubiq K"), ("glyc", "N-Glyc N")]
# pure-positive share of the threshold-sampled test partitions (analysis/crosseval_summary.py, PP)
PP = {"phosphorylation_y": 3.5, "phosphorylation_st": 3.8, "ubiquitination_k": 12.6,
      "sumoylation_k": 16.5, "acetylation_k": 20.3, "methylation_k": 27.8,
      "methylation_r": 31.6, "glycosylation_n": 35.9}


def exact_p(x, y):
    """Two-sided exact permutation p for Spearman's rho (all n! pairings)."""
    rx = rankdata(x); ry = np.array(list(permutations(rankdata(y))))
    rx = rx - rx.mean(); ry = ry - ry.mean(axis=1, keepdims=True)
    r = ry @ rx / np.sqrt((ry ** 2).sum(axis=1) * (rx ** 2).sum())
    return float(np.mean(np.abs(r) >= abs(spearmanr(x, y).statistic) - 1e-9))


def label_rate(ax):
    d = pd.read_csv(S3 / "depth_plot_splits.tsv", sep="\t")
    d = d[(d.family == "esm") & (d.train == "rebuilt")]
    for t in TASKS:
        g = d[d.task == t].groupby("depth_bin")[["depth_median", "pos_rate"]].mean()
        if t == "acetylation_k":
            ax.plot(g.depth_median, 100 * g.pos_rate, "-o", color=INK, ms=3.5, lw=1.4, zorder=4)
            ax.text(g.depth_median.iloc[-1] * 1.15, 100 * g.pos_rate.iloc[-1], "Acetyl K", fontsize=6, va="center")
        elif t == "phosphorylation_y":
            ax.plot(g.depth_median, 100 * g.pos_rate, "--o", color=INK2, ms=3, lw=1.1, mfc="white", zorder=3)
            ax.text(g.depth_median.iloc[-1] * 1.15, 100 * g.pos_rate.iloc[-1], "Phospho Y", fontsize=6,
                    va="center", color=INK2)
        else:
            ax.plot(g.depth_median, 100 * g.pos_rate, "-", color=AXIS, lw=1.0, zorder=2)
    ax.set_xscale("log")
    ax.set_xlabel("Merged annotation depth (bin median)")
    ax.set_ylabel("Positives among retained sites,\npooled within bin (%)")
    ax.set_ylim(0, 40)
    ax.set_title("Unrestricted reconstruction", loc="left", fontsize=6.5, color=INK2)
    hgrid(ax)


def probe_scatter(ax):
    o = pd.read_csv(S3 / "ridge_acetylation_oof.tsv", sep="\t")
    r = pd.read_csv(S3 / "ridge_summary.tsv", sep="\t").set_index(["group", "family"]).raw_r2
    for fam in ("ppi", "esm"):
        g = o[o.family == fam]
        ax.scatter(g.observed, g.predicted, s=1.5, alpha=0.18, color=CH[FAM[fam]], linewidth=0,
                   rasterized=True)
    lo, hi = -0.1, max(o.observed.max(), o.predicted.max()) + 0.1
    ax.plot([lo, hi], [lo, hi], color=INK2, lw=0.7, ls=(0, (2, 2)))
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("Observed log(1 + depth)")
    ax.set_ylabel("Held-out prediction")
    h = [Line2D([], [], ls="", marker="o", color=CH[FAM[f]], ms=4,
                label=f"{CH_LABEL[FAM[f]]}  R² = {r[('acet', f)]:.2f}") for f in ("ppi", "esm")]
    ax.legend(handles=h, loc="upper left", handletextpad=0.2, borderaxespad=0.1)
    ax.set_title("Acetylation K, 6,249 proteins", loc="left", fontsize=6.5, color=INK2)


def probe_r2(ax):
    r = pd.read_csv(S3 / "ridge_summary.tsv", sep="\t").set_index(["group", "family"])
    y = np.arange(len(GROUP))[::-1]
    for k, fam in enumerate(("ppi", "esm")):
        yy = y + (0.15 if k == 0 else -0.15)
        raw = [r.loc[(g, fam), "raw_r2"] for g, _ in GROUP]
        res = [r.loc[(g, fam), "residual_r2"] for g, _ in GROUP]
        for a, b, yv in zip(raw, res, yy):
            ax.plot([b, a], [yv, yv], color=CH[FAM[fam]], lw=1.0, alpha=0.6, zorder=2)
        ax.scatter(raw, yy, s=20, color=CH[FAM[fam]], edgecolor="white", linewidth=0.6, zorder=3)
        ax.scatter(res, yy, s=20, facecolor="white", edgecolor=CH[FAM[fam]], linewidth=1.0, zorder=3)
    ax.set_yticks(y, [n for _, n in GROUP])
    ax.set_xlim(0, 0.62)
    ax.set_xlabel("Held-out R²")
    h = [Line2D([], [], ls="", marker="o", color=INK2, ms=4, label="raw target: log(1 + depth)"),
         Line2D([], [], ls="", marker="o", mfc="white", mec=INK2, ms=4,
                label="residual target: after candidate-count adjustment")]
    ax.legend(handles=h, loc="lower left", bbox_to_anchor=(-0.05, 1.0), ncol=1, handletextpad=0.2,
              borderaxespad=0.1)
    vgrid(ax)


def dumbbell(ax, before, after, title_before, xlabel, ref_color=INK2):
    """Per task, per feature family: before -> after."""
    y = np.arange(len(TASKS))[::-1]
    ax.axvline(0, color=INK2, lw=0.7)
    for k, fam in enumerate(("ppi", "esm")):
        yy = y + (0.17 if k == 0 else -0.17)
        b = [before[(t, fam)] for t in TASKS]
        a = [after[(t, fam)] for t in TASKS]
        for bv, av, yv in zip(b, a, yy):
            ax.annotate("", xy=(av, yv), xytext=(bv, yv),
                        arrowprops=dict(arrowstyle="-|>", color=CH[FAM[fam]], lw=0.9, shrinkA=2.5,
                                        shrinkB=3, mutation_scale=6, alpha=0.8))
        ax.scatter(b, yy, s=14, facecolor="white", edgecolor=ref_color, linewidth=0.9, zorder=3)
        ax.scatter(a, yy, s=16, color=CH[FAM[fam]], edgecolor="white", linewidth=0.5, zorder=3)
    ax.set_yticks(y, [SHORT[t] for t in TASKS])
    ax.set_xlabel(xlabel)
    vgrid(ax)


def score_shift(ax):
    s = pd.read_csv(S3 / "score_depth_summary.tsv", sep="\t")
    before = {(r.task, r.family): r.baseline for r in s.itertuples()}
    after = {(r.task, r.family): r.augmented for r in s.itertuples()}
    dumbbell(ax, before, after, "sequence only", "Spearman ρ, mean site score vs depth")
    ax.set_title("Threshold training, unrestricted test", loc="left", fontsize=6.5, color=INK2)
    h = [Line2D([], [], ls="", marker="o", mfc="white", mec=INK2, ms=4, label="sequence only"),
         Line2D([], [], ls="", marker="o", color=CH["interaction"], ms=4, label="+ interaction"),
         Line2D([], [], ls="", marker="o", color=CH["esm2"], ms=4, label="+ ESM-2")]
    ax.legend(handles=h, loc="upper left", bbox_to_anchor=(-0.05, -0.2), ncol=1, handletextpad=0.1,
              columnspacing=0.6, borderaxespad=0)


def neighbours(ax):
    k = pd.read_csv(S3 / "knn_summary.tsv", sep="\t")
    k = k[k.family == "ppi"].set_index("task")
    y = np.arange(len(TASKS))[::-1]
    ax.axvline(0, color=INK2, lw=0.7)
    for t, yv in zip(TASKS, y):
        a, b = k.loc[t, "augmented_replica"], k.loc[t, "augmented_rebuilt"]
        ax.plot([a, b], [yv, yv], color=AXIS, lw=1.2, zorder=2)
        ax.scatter(a, yv, s=20, color=CH["interaction"], edgecolor="white", linewidth=0.6, zorder=3)
        ax.scatter(b, yv, s=20, facecolor="white", edgecolor=CH["interaction"], linewidth=1.0, zorder=3)
    ax.set_yticks(y, [SHORT[t] for t in TASKS])
    ax.set_xlabel("ρ, score vs positive fraction of 25 neighbours")
    ax.set_title("Interaction space, same neighbours", loc="left", fontsize=6.5, color=INK2)
    h = [Line2D([], [], ls="", marker="o", color=CH["interaction"], ms=4, label="threshold-sampled labels"),
         Line2D([], [], ls="", marker="o", mfc="white", mec=CH["interaction"], ms=4, label="unrestricted labels")]
    ax.legend(handles=h, loc="upper left", bbox_to_anchor=(-0.05, -0.2), ncol=1, handletextpad=0.1,
              columnspacing=0.6, borderaxespad=0)
    vgrid(ax)


def depth_bins(ax):
    d = pd.read_csv(S3 / "depth_plot_splits.tsv", sep="\t")
    d = d[(d.family == "esm") & (d.task == "acetylation_k")]
    for tr, style, mfc, lab in [("replica", "-", INK, "threshold training"),
                                ("rebuilt", "--", "white", "unrestricted training")]:
        g = d[d.train == tr].groupby("depth_bin")[["depth_median", "d_pct", "pos_rate"]].mean()
        ax.plot(g.depth_median, 100 * g.d_pct, style, marker="o", color=INK, mfc=mfc, ms=3.5, lw=1.2, label=lab)
    ax.axhline(0, color=INK2, lw=0.7)
    ax.set_xscale("log")
    ax.set_xlabel("Merged annotation depth (bin median)")
    ax.set_ylabel("Change in score percentile (points)")
    ax.set_title("Acetylation K, + ESM-2, unrestricted test", loc="left", fontsize=6.5, color=INK2)
    ax.legend(loc="lower left", handlelength=1.8, borderaxespad=0.1)
    hgrid(ax)


def dose(axes):
    v = pd.read_csv(RES / "crosseval_verify.tsv", sep="\t")
    v = v[v.family == "interaction"]
    cells = [(("replica", "replica"), "TT"), (("replica", "rebuilt"), "TU"), (("rebuilt", "rebuilt"), "UU")]
    keep = [t for t in TASKS if t != "glycosylation_n"]
    for ax, ((tr, ev), name) in zip(axes, cells):
        m = v[(v.train == tr) & (v["eval"] == ev)].groupby("task").d_auroc.mean()
        x7, y7 = [PP[t] for t in keep], [m[t] for t in keep]
        rho, p = spearmanr(x7, y7).statistic, exact_p(x7, y7)
        ax.axhline(0, color=INK2, lw=0.7)
        ax.scatter(x7, y7, s=20, color=CH["interaction"], edgecolor="white", linewidth=0.6, zorder=3)
        ax.scatter(PP["glycosylation_n"], m["glycosylation_n"], s=20, facecolor="white",
                   edgecolor=CH["interaction"], linewidth=1.0, zorder=3)
        ax.set_title(f"{name}   ρ = {rho:+.2f}, p = {p:.3f} (n = 7)", loc="left", fontsize=6.5)
        ax.set_xlabel("Pure-positive share of the\nthreshold-sampled test partition (%)")
        hgrid(ax)
    axes[0].set_ylabel("ΔAUROC, interaction")
    axes[2].text(0.98, 0.04, "open: N-Glyc N, excluded from ρ", transform=axes[2].transAxes,
                 ha="right", fontsize=5.6, color=INK2)


def main():
    setup()
    fig = plt.figure(figsize=(180 * MM, 215 * MM))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1.05, 0.8], hspace=0.85, wspace=0.6)
    a = fig.add_subplot(gs[0, 0]); label_rate(a); label(a, "a", x=-0.3, y=1.06)
    b = fig.add_subplot(gs[0, 1]); probe_scatter(b); label(b, "b", x=-0.3, y=1.06)
    c = fig.add_subplot(gs[0, 2]); probe_r2(c); label(c, "c", x=-0.55, y=1.06)
    d = fig.add_subplot(gs[1, 0]); score_shift(d); label(d, "d", x=-0.42, y=1.06)
    e = fig.add_subplot(gs[1, 1]); neighbours(e); label(e, "e", x=-0.42, y=1.06)
    f = fig.add_subplot(gs[1, 2]); depth_bins(f); label(f, "f", x=-0.3, y=1.06)
    sub = gs[2, :].subgridspec(1, 3, wspace=0.12)
    g_axes = [fig.add_subplot(sub[0, i]) for i in range(3)]
    dose(g_axes)
    lo = min(ax.get_ylim()[0] for ax in g_axes) - 0.01
    hi = max(ax.get_ylim()[1] for ax in g_axes) + 0.02
    for i, ax in enumerate(g_axes):
        ax.set_ylim(lo, hi)
        if i:
            ax.tick_params(labelleft=False)
    label(g_axes[0], "g", x=-0.3)
    save(fig, "Figure3_explanation")


if __name__ == "__main__":
    main()
