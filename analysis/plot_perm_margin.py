#!/usr/bin/env python3
"""Scatter of the permutation margin against test-partition homogeneity.

The main text reports a rank correlation that is -0.43 over eight tasks and
-0.89 once N-glycosylation is excluded. Stating that in prose invites the
reading that a point was dropped to obtain the second number. The scatter shows
why it is set aside: N-glycosylation sits at the edge of the baseline axis, its
sequon-only classifier already reaching 0.9428, so there is almost no room for a
protein-level channel to move it in either direction. The reader can see the
geometry instead of taking the exclusion on trust.

    python plot_perm_margin.py -o figure_S5.png
"""
import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

# real minus permuted, interaction embedding, threshold-trained, natural
# negatives; means over three partitions. From crosseval_verify.tsv.
MARGIN = {"Phospho S/T": +0.0026, "Phospho Y": +0.0208, "Acetyl K": -0.1285,
          "Meth K/R": -0.0613, "Meth R": -0.0706, "Sumo K": -0.0558,
          "Ubiq K": -0.0033, "N-Glyc N": +0.0074}
# pure-positive share of the rebuilt-evaluated test partition, per cent
PUREPOS = {"Phospho Y": 3.5, "Phospho S/T": 3.8, "Ubiq K": 12.6, "Sumo K": 16.5,
           "Acetyl K": 20.3, "Meth K/R": 27.8, "Meth R": 31.6, "N-Glyc N": 35.9}
# sequence-only AUROC on the protein-disjoint reconstruction, natural negatives
BASELINE = {"Phospho S/T": 0.8712, "Phospho Y": 0.8005, "Acetyl K": 0.7995,
            "Meth K/R": 0.7356, "Meth R": 0.7678, "Sumo K": 0.7710,
            "Ubiq K": 0.7851, "N-Glyc N": 0.9615}
EXCLUDED = "N-Glyc N"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default="figure_S5.png")
    ap.add_argument("--dpi", type=int, default=300)
    a = ap.parse_args()

    tasks = list(MARGIN)
    keep = [t for t in tasks if t != EXCLUDED]
    r8 = spearmanr([PUREPOS[t] for t in tasks], [MARGIN[t] for t in tasks])
    r7 = spearmanr([PUREPOS[t] for t in keep], [MARGIN[t] for t in keep])

    fig, ax = plt.subplots(1, 2, figsize=(6.8, 3.0))

    # Labels are placed away from the zero line, which several points sit on:
    # a label at a fixed offset lands on the dashed rule and is unreadable.
    def place(a, x, y, t, xs):
        # Above the marker for a positive margin, below for a negative one, so no
        # label lands on the zero rule that several points sit against; and to the
        # left for points near the right edge, which would otherwise run off it.
        lo, hi = min(xs), max(xs)
        right = x < lo + 0.72 * (hi - lo)
        a.annotate(t, (x, y), fontsize=5.6,
                   xytext=(4 if right else -4, 5 if y >= 0 else -9),
                   ha="left" if right else "right",
                   va="bottom" if y >= 0 else "top",
                   textcoords="offset points")

    xs0 = list(PUREPOS.values())
    for t in tasks:
        ex = t == EXCLUDED
        ax[0].scatter(PUREPOS[t], MARGIN[t], s=34,
                      facecolor="white" if ex else "#3b6ea5",
                      edgecolor="#c0392b" if ex else "#3b6ea5",
                      linewidth=1.4 if ex else 0, zorder=3)
        place(ax[0], PUREPOS[t], MARGIN[t], t, xs0)
    ax[0].axhline(0, color="k", lw=.6, ls="--")
    ax[0].set_xlabel("pure-positive share of the test partition (%)", fontsize=8)
    ax[0].set_ylabel("real minus permuted, $\\Delta$AUROC", fontsize=8)
    fmt = lambda v: f"{v:+.2f}".replace("-", "\u2212")
    ax[0].set_title(f"a  $\\rho$ = {fmt(r8.statistic)} over eight tasks, "
                    f"{fmt(r7.statistic)} over seven", loc="left", fontsize=8)

    xs1 = list(BASELINE.values())
    for t in tasks:
        ex = t == EXCLUDED
        ax[1].scatter(BASELINE[t], MARGIN[t], s=34,
                      facecolor="white" if ex else "#3b6ea5",
                      edgecolor="#c0392b" if ex else "#3b6ea5",
                      linewidth=1.4 if ex else 0, zorder=3)
        place(ax[1], BASELINE[t], MARGIN[t], t, xs1)
    ax[1].axhline(0, color="k", lw=.6, ls="--")
    ax[1].set_xlabel("sequence-only AUROC", fontsize=8)
    ax[1].set_title("b  why N-glycosylation is set aside", loc="left", fontsize=8)

    for x in ax:
        x.margins(y=0.16, x=0.10)
        x.tick_params(labelsize=7)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(a.out, dpi=a.dpi, bbox_inches="tight")
    print(f"wrote {a.out}")
    print(f"  rho over eight tasks {r8.statistic:+.3f}  p = {r8.pvalue:.4f}")
    print(f"  rho excluding {EXCLUDED} {r7.statistic:+.3f}  p = {r7.pvalue:.4f}")


if __name__ == "__main__":
    sys.exit(main())
