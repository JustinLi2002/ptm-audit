"""Figure 4. The reversal in a synthetic system, and sensitivity to missing positives."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figstyle import CH, INK, INK2, MM, RAMP5, REPO, S3, hgrid, label, save, setup

SHOW_T = [1, 2, 5, 10, 20]          # on the ordinal ramp; T = 0 (no threshold) is drawn as the reference
MERGE = {"methylation_k": "methyl", "methylation_r": "methyl"}


def sweep():
    d = pd.read_csv(REPO / "figures" / "figure6_sweep.tsv", sep="\t")
    return d.groupby(["threshold", "r2_target"])[["r2_measured", "d_real_natural", "rho_delta"]].mean().reset_index()


def synth_panel(ax, g, col, ylabel):
    ax.axhline(0, color=INK2, lw=0.7)
    s = g[g.threshold == 0].sort_values("r2_measured")
    ax.plot(s.r2_measured, s[col], "--o", color=INK2, ms=3, lw=1.2, mfc="white", label="T = 0 (no threshold)")
    for c, t in zip(RAMP5, SHOW_T):
        s = g[g.threshold == t].sort_values("r2_measured")
        ax.plot(s.r2_measured, s[col], "-o", color=c, ms=3.2, lw=1.4, mec="white", mew=0.5, label=f"T = {t}")
    ax.set_xlabel("Feature fidelity\n(ridge R² from feature to log depth)")
    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.01, g.r2_measured.max() + 0.02)
    hgrid(ax)


def fn_panel(ax):
    ax.axhline(0, color=INK2, lw=0.7)
    for arm, fam, name in [("ppi", "interaction", "Interaction"), ("esm", "esm2", "ESM-2")]:
        d = pd.read_csv(S3 / f"fn_sensitivity_{arm}.tsv", sep="\t")
        d = d[d["alpha"].notna() & d["target_inflation"].isna()].copy()
        d["unit"] = d["task"].map(MERGE).fillna(d["task"])
        curve = d.groupby(["alpha", "unit"]).d_auroc.mean().groupby("alpha").mean()
        ax.plot(curve.index, curve.values, "-o", color=CH[fam], ms=3.5, lw=1.4, mec="white", mew=0.5)
        a0, a1 = curve.loc[0.0], curve.loc[1.0]
        ax.text(1.03, a1, name, fontsize=6, va="center", color=INK)
        up = arm == "esm"
        ax.text(0.04, a0 + (0.004 if up else -0.004), f"{a0 / a1:.0%} remains",
                fontsize=5.8, color=INK2, va="bottom" if up else "top")
        print(arm, curve.round(4).to_dict())
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xlim(-0.03, 1.2)
    ax.set_xlabel("α  (0: depth–rate relationship all study bias;\n1: annotation complete)")
    ax.set_ylabel("Mean ΔAUROC, unrestricted test (7 units)")
    ax.set_title("Threshold training; sensitivity model", loc="left", fontsize=6.5, color=INK2)
    hgrid(ax)


def main():
    setup()
    g = sweep()
    fig, axes = plt.subplots(1, 3, figsize=(180 * MM, 62 * MM), gridspec_kw={"wspace": 0.42})
    synth_panel(axes[0], g, "d_real_natural", "ΔAUROC from the entity-constant channel")
    axes[0].set_title("Synthetic: natural test", loc="left", fontsize=6.5, color=INK2)
    synth_panel(axes[1], g, "rho_delta", "Change in ρ(score, depth)")
    axes[1].set_title("Synthetic: score–depth relationship", loc="left", fontsize=6.5, color=INK2)
    hs, ls = axes[0].get_legend_handles_labels()
    fig.legend(hs, ls, loc="upper center", bbox_to_anchor=(0.36, -0.13), ncol=6, fontsize=6,
               handlelength=1.8, columnspacing=1.2)
    fn_panel(axes[2])
    for ax, l in zip(axes, "abc"):
        label(ax, l, x=-0.3, y=1.06)
    save(fig, "Figure4_synthetic_sensitivity")
    print(g[g.threshold.isin([5, 10, 20])].round(3).to_string())


if __name__ == "__main__":
    main()
