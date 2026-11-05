"""Supplementary Fig. 3. Dose-response in all four cells, for all three protein-level features."""
from itertools import permutations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from figstyle import CH, CH_LABEL, INK2, MM, RES, hgrid, label, save, setup

# pure-positive share of the threshold-sampled test partitions (analysis/crosseval_summary.py, PP)
PP = {"phosphorylation_y": 3.5, "phosphorylation_st": 3.8, "ubiquitination_k": 12.6,
      "sumoylation_k": 16.5, "acetylation_k": 20.3, "methylation_k": 27.8,
      "methylation_r": 31.6, "glycosylation_n": 35.9}
CELLS = [("replica", "replica", "TT"), ("replica", "rebuilt", "TU"),
         ("rebuilt", "replica", "UT"), ("rebuilt", "rebuilt", "UU")]
FAMS = ["interaction", "esm2", "prott5"]


def exact_p(x, y):
    """Two-sided exact permutation p for Spearman's rho (all n! pairings)."""
    rx = rankdata(x); ry = np.array(list(permutations(rankdata(y))))
    rx = rx - rx.mean(); ry = ry - ry.mean(axis=1, keepdims=True)
    r = ry @ rx / np.sqrt((ry ** 2).sum(axis=1) * (rx ** 2).sum())
    return float(np.mean(np.abs(r) >= abs(spearmanr(x, y).statistic) - 1e-9))


def fmt(r, p):
    return f"ρ = {r:+.2f}, p = {p:.3f}".replace("-", "−").replace("+0.00", "0.00")


def main():
    setup()
    v = pd.read_csv(RES / "crosseval_verify.tsv", sep="\t")
    keep = [t for t in PP if t != "glycosylation_n"]
    fig, axes = plt.subplots(3, 4, figsize=(180 * MM, 135 * MM), sharex=True, sharey=True,
                             gridspec_kw={"hspace": 0.62, "wspace": 0.12})
    for i, fam in enumerate(FAMS):
        for j, (tr, ev, name) in enumerate(CELLS):
            ax = axes[i, j]
            m = v[(v.family == fam) & (v.train == tr) & (v["eval"] == ev)].groupby("task").d_auroc.mean()
            x7, y7 = [PP[t] for t in keep], [m[t] for t in keep]
            x8, y8 = list(PP.values()), [m[t] for t in PP]
            r7, p7 = spearmanr(x7, y7).statistic, exact_p(x7, y7)
            r8, p8 = spearmanr(x8, y8).statistic, exact_p(x8, y8)
            ax.axhline(0, color=INK2, lw=0.7)
            ax.scatter(x7, y7, s=18, color=CH[fam], edgecolor="white", linewidth=0.6, zorder=3)
            ax.scatter(PP["glycosylation_n"], m["glycosylation_n"], s=18, facecolor="white",
                       edgecolor=CH[fam], linewidth=1.0, zorder=3)
            ax.set_title(f"7 tasks: {fmt(r7, p7)}", loc="left", fontsize=6, pad=9)
            ax.text(0, 1.02, f"8 tasks: {fmt(r8, p8)}", transform=ax.transAxes,
                    ha="left", va="bottom", fontsize=5.6, color=INK2)
            if i == 0:
                ax.text(0.5, 1.3, name, transform=ax.transAxes, ha="center", va="bottom",
                        fontsize=8, fontweight="bold")
            hgrid(ax)
            print(f"{fam:12s} {name}  n7 {r7:+.3f} p {p7:.4f} | n8 {r8:+.3f} p {p8:.4f}")
        axes[i, 0].set_ylabel(f"ΔAUROC, {CH_LABEL[fam]}")
        label(axes[i, 0], "abc"[i], x=-0.42, y=1.08)
    fig.supxlabel("Pure-positive share of the threshold-sampled test partition (%)", fontsize=7, y=0.045)
    save(fig, "SuppFig3_dose_response")


if __name__ == "__main__":
    main()
