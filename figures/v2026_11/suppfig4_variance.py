"""Supplementary Fig. 4. Within-protein variance ratio and ICC in all four cells, for all three features,
with the earlier Fig. 4a dot plot as a worked example."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import pandas as pd

from figstyle import CH, CH_LABEL, INK, INK2, MM, MUTED, OUT, RES, hgrid, label, save, setup

# pure-positive share of the threshold-sampled test partitions (analysis/crosseval_summary.py, PP)
PP = {"phosphorylation_y": 3.5, "phosphorylation_st": 3.8, "ubiquitination_k": 12.6,
      "sumoylation_k": 16.5, "acetylation_k": 20.3, "methylation_k": 27.8,
      "methylation_r": 31.6, "glycosylation_n": 35.9}
FAMS = ["interaction", "esm2", "prott5"]
# training construction sets the colour (feature colour vs grey); test construction sets line and fill
CELLS = [("replica", "replica", "TT", "--", False), ("replica", "rebuilt", "TU", "-", True),
         ("rebuilt", "replica", "UT", "--", False), ("rebuilt", "rebuilt", "UU", "-", True)]
# per-site predictions for the worked example, fetched from HRP/pdisjoint_runs_v2 on 2026-11-01
PRED = OUT / "data_hpc_20261101" / "acetylation_k__replica__{}__split0__on_rebuilt.pred.tsv"


def load():
    d = pd.read_csv(RES / "icc_audit.tsv", sep="\t").drop_duplicates(["task", "train", "eval", "split", "family"])
    # ratio of partition-mean mean squares (as in results/icc_*.txt); ICC is the partition mean
    g = d.groupby(["family", "train", "eval", "task"])[["msw_base", "msw", "icc"]].mean()
    g["ratio"] = g.msw_base / g.msw
    return g.reset_index()


def msw(d):
    dev = d.y_pred - d.groupby("protein").y_pred.transform("mean")
    return (dev ** 2).sum() / (len(d) - d.protein.nunique())


def example(ax):
    """Old Fig. 4a: acetylation K, TU, interaction feature, first partition, 18 proteins spread over the range."""
    fr = {c: pd.read_csv(str(PRED).format(c), sep="\t") for c in ("baseline", "ppi")}
    k = fr["baseline"].groupby("protein").size()
    prots = k[k >= 5].index
    fr = {c: d[d.protein.isin(prots)] for c, d in fr.items()}
    ranked = fr["ppi"].groupby("protein").y_pred.mean().sort_values().index
    show = list(ranked[::max(1, len(ranked) // 18)])[:18]
    rng = np.random.default_rng(0)
    for off, cond, col in [(0, "baseline", MUTED), (19, "ppi", CH["interaction"])]:
        d = fr[cond]
        for i, p in enumerate(show):
            v = d[d.protein == p].y_pred.values
            ax.scatter(i + off + rng.uniform(-0.2, 0.2, len(v)), v, s=2.5, alpha=0.55, color=col,
                       linewidths=0, rasterized=True)
            ax.plot([i + off - 0.38, i + off + 0.38], [v.mean()] * 2, color=INK, lw=0.8)
    r = msw(fr["baseline"]) / msw(fr["ppi"])
    print(f"example: {len(prots)} proteins, MSW ratio {r:.2f}")
    ax.set_xticks([8.5, 27.5], ["Sequence only", "With interaction feature"])
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-1, 37)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Predicted score")
    ax.set_title(f"Acetylation K, TU, first partition: 18 of {len(prots)} test proteins shown "
                 f"(variance ratio {r:.2f})", loc="left", fontsize=7)
    hgrid(ax)


def main():
    setup()
    g = load()
    order = sorted(PP, key=PP.get)
    fig = plt.figure(figsize=(180 * MM, 160 * MM))
    gs = fig.add_gridspec(3, 3, height_ratios=[0.8, 1, 1], hspace=0.42, wspace=0.28)
    ax0 = fig.add_subplot(gs[0, :])
    example(ax0)
    axes = np.array([[fig.add_subplot(gs[r, c]) for c in range(3)] for r in (1, 2)])
    for j, fam in enumerate(FAMS):
        for tr, ev, name, ls, filled in CELLS:
            col = CH[fam] if tr == "replica" else MUTED
            s = g[(g.family == fam) & (g.train == tr) & (g["eval"] == ev)].set_index("task").loc[order]
            x = [PP[t] for t in order]
            for row, key in enumerate(("ratio", "icc")):
                axes[row, j].plot(x, s[key], ls=ls, color=col, lw=1.1, marker="o", ms=3.4,
                                  mfc=col if filled else "white", mec=col, mew=0.9, zorder=3)
            print(f"{fam:12s} {name} ratio {s.ratio.min():.2f}-{s.ratio.max():.2f} "
                  f"icc {s.icc.min():.2f}-{s.icc.max():.2f}")
        axes[0, j].axhline(1, color=INK2, lw=0.7)
        axes[0, j].set_yscale("log")
        axes[0, j].set_yticks([0.8, 1, 2, 4, 8], ["0.8", "1", "2", "4", "8"])
        axes[0, j].minorticks_off()
        axes[0, j].set_title(CH_LABEL[fam], loc="left", fontsize=7.5)
        axes[0, j].sharex(axes[1, j])
        axes[0, j].tick_params(labelbottom=False)
        axes[1, j].set_ylim(0, 1)
        axes[1, j].set_xlabel("Pure-positive share of the\nthreshold-sampled test partition (%)")
        for row in range(2):
            hgrid(axes[row, j])
    axes[0, 0].set_ylabel("Within-protein mean square,\nsequence-only / augmented")
    axes[1, 0].set_ylabel("ICC(1,1), augmented model")
    label(ax0, "a", x=-0.098, y=1.06)
    label(axes[0, 0], "b", x=-0.36, y=1.06)
    label(axes[1, 0], "c", x=-0.36, y=1.06)
    hl = [Line2D([], [], ls=ls, color=INK if tr == "replica" else MUTED, marker="o", ms=3.4,
                 mfc=(INK if tr == "replica" else MUTED) if filled else "white",
                 mec=INK if tr == "replica" else MUTED, label=name)
          for tr, ev, name, ls, filled in CELLS]
    fig.legend(handles=hl, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=4, handlelength=2.6,
               columnspacing=1.6)
    fig.subplots_adjust(bottom=0.11)
    save(fig, "SuppFig4_variance")


if __name__ == "__main__":
    main()
