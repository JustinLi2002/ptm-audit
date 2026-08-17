#!/usr/bin/env python3
"""Figure for the depth-stratified dose-response.

    python plot_depth_strata.py depth_esm2_thr.tsv depth_esm2_unr.tsv \
        -o figure_S4.png

Two rows of panels, one task per column: the change in mean score percentile by
annotation-depth bin and candidate-count bin, under threshold-sampled training
(top) and unrestricted training (bottom, the negative control). A depth effect
is a vertical gradient that persists in every column of a panel; a protein-size
effect is a horizontal one. Both are present, and the note in the caption
explains why: the donor threshold is defined on a site count, which grows with
protein size, so the two gradients are projections of one rule.
"""
import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

SHORT = {"phosphorylation_st": "Phospho S/T", "phosphorylation_y": "Phospho Y",
         "acetylation_k": "Acetyl K", "methylation_k": "Meth K/R",
         "methylation_r": "Meth R", "sumoylation_k": "Sumo K",
         "ubiquitination_k": "Ubiq K", "glycosylation_n": "N-Glyc N"}
ORDER = ["phosphorylation_st", "phosphorylation_y", "acetylation_k",
         "methylation_k", "methylation_r", "sumoylation_k",
         "ubiquitination_k", "glycosylation_n"]


def grid(df, task):
    s = df[df.task == task]
    return s.pivot(index="depth_bin", columns="size_bin", values="d_pct")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("threshold_tsv")
    ap.add_argument("unrestricted_tsv")
    ap.add_argument("-o", "--out", default="figure_depth_strata.png")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    thr = pd.read_csv(args.threshold_tsv, sep="\t")
    unr = pd.read_csv(args.unrestricted_tsv, sep="\t")
    tasks = [t for t in ORDER if t in set(thr.task)]

    vmax = max(np.nanmax(np.abs(thr.d_pct)), np.nanmax(np.abs(unr.d_pct)))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    fig, axes = plt.subplots(2, len(tasks), figsize=(1.55 * len(tasks), 4.4),
                             squeeze=False)
    for j, task in enumerate(tasks):
        for i, df in enumerate((thr, unr)):
            ax = axes[i][j]
            g = grid(df, task)
            im = ax.imshow(g.to_numpy(), cmap="RdBu_r", norm=norm,
                           aspect="auto", origin="upper")
            ax.set_xticks([])
            # every panel carries its own tick marks: depth is binned by
            # quantile and several tasks admit only three or four bins because
            # depth is heavily tied at one and two sites. A shared 0-4 axis
            # would imply five bins everywhere.
            ax.set_yticks(range(len(g)))
            ax.set_yticklabels([str(k + 1) for k in g.index], fontsize=6)
            ax.tick_params(length=1.5, pad=1)
            if i == 0:
                ax.set_title(f"{SHORT.get(task, task)}\n({len(g)} depth bins)",
                             fontsize=7, pad=4)
            if j == 0:
                ax.set_ylabel("threshold\ntraining" if i == 0
                              else "unrestricted\ntraining", fontsize=7)
    for ax in axes[1]:
        ax.set_xlabel("candidate\ncount →", fontsize=6)
    axes[0][0].text(-0.75, 0.5, "annotation depth →", rotation=90, fontsize=6,
                    va="center", ha="center", transform=axes[0][0].transAxes)

    cb = fig.colorbar(im, ax=axes, fraction=0.015, pad=0.012)
    cb.set_label("Δ mean score percentile", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    print(f"wrote {args.out}")

    print("\nspan of the depth gradient, pooled over sizes "
          "(shallowest bin minus deepest):")
    for task in tasks:
        f = lambda d: (d[(d.task == task)].groupby("depth_bin")
                       .apply(lambda g: np.average(g.d_pct.fillna(0),
                                                   weights=g.n_sites)))
        a, b = f(thr), f(unr)
        print(f"  {SHORT.get(task, task):13s} threshold {a.iloc[0] - a.iloc[-1]:+.3f}"
              f"   unrestricted {b.iloc[0] - b.iloc[-1]:+.3f}")


if __name__ == "__main__":
    sys.exit(main())
