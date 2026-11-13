"""Figure 4 (v5). The reversal in a synthetic system, and a ratio-matched control.
v5 differs from v4 only in panel c: the tick label reads "Ratio-matched" (was "Matched"), as in Fig. 2g,h,
and panel c is 15% wider, with less margin beside the outer groups, so that the tick labels do not overlap.
Panels a and b are those of fig4_synthetic_sensitivity.py; the missing-positive panel moved to Supplementary Fig. 5."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fig4_synthetic_sensitivity import sweep, synth_panel
from figstyle import INK, INK2, MM, MUTED, RAMP5, hgrid, label, save, setup

RUNS = Path(__file__).resolve().parents[2] / "ratio_matched" / "synthetic" / "runs.tsv"
GROUPS = [("T", "Threshold", RAMP5[3], "o", RAMP5[3]),
          ("U_matched", "Ratio-matched", INK2, "s", INK2),
          ("U", "Unrestricted", INK2, "o", "white")]


def matched_panel(ax):
    r = pd.read_csv(RUNS, sep="\t")
    p = r.pivot(index="seed", columns="construction", values="dauroc_real_natural")
    rate = r.groupby("construction").pos_rate.mean()
    x = np.arange(len(GROUPS))
    ax.axhline(0, color=INK2, lw=0.7)
    for seed in p.index:
        ax.plot(x, [p.loc[seed, g] for g, *_ in GROUPS], color=MUTED, lw=0.6, alpha=0.7, zorder=1)
    for xi, (g, _, c, m, fc) in zip(x, GROUPS):
        ax.scatter(np.full(len(p), xi), p[g], s=16, marker=m, color=c, facecolor=fc, linewidth=0.9, zorder=3)
        ax.plot([xi - 0.22, xi + 0.22], [p[g].mean()] * 2, color=INK, lw=1.4, zorder=4)
    ax.set_xticks(x, [f"{lab}\n{100 * rate[g]:.0f}%" for g, lab, *_ in GROUPS], fontsize=6.5)
    ax.set_xlabel("Training set (positive rate)")
    ax.set_xlim(-0.38, len(GROUPS) - 0.62)
    ax.set_ylabel("ΔAUROC from the entity-constant feature")
    meas = r.r2_measured.mean()
    ax.set_title(f"Ratio-matched, T = 10, R² = {meas:.2f}", loc="left", fontsize=6.5, color=INK2)
    hgrid(ax)
    print(p.round(4).to_string())


def main():
    setup()
    g = sweep()
    fig, axes = plt.subplots(1, 3, figsize=(180 * MM, 62 * MM), gridspec_kw={"wspace": 0.42, "width_ratios": [1, 1, 1.15]})
    synth_panel(axes[0], g, "d_real_natural", "ΔAUROC from the entity-constant feature")
    axes[0].set_title("Synthetic: test without threshold", loc="left", fontsize=6.5, color=INK2)
    synth_panel(axes[1], g, "rho_delta", "Change in ρ(score, depth)")
    axes[1].set_title("Synthetic: score–depth relationship", loc="left", fontsize=6.5, color=INK2)
    hs, ls = axes[0].get_legend_handles_labels()
    fig.legend(hs, ls, loc="upper center", bbox_to_anchor=(0.36, -0.13), ncol=6, fontsize=6,
               handlelength=1.8, columnspacing=1.2)
    matched_panel(axes[2])
    for ax, l in zip(axes, "abc"):
        label(ax, l, x=-0.3, y=1.06)
    save(fig, "Figure4_synthetic_matched_v5")


if __name__ == "__main__":
    main()
