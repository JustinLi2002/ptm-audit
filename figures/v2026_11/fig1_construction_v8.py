"""Figure 1 (v8). Panels a-d exactly as v7; a new full-width panel e, a schematic of what each construction
associates with annotation depth (after Fig. 1d of the 2026-09-21 draft, relabelled: annotated fraction, not
modification rate; the threshold side drawn as the step the rule produces). Shapes follow the data of Fig. 3a
(acetylation K: every retained site positive below depth ten, then 0.19-0.37; unrestricted 0.04 rising to 0.48),
but the panel carries no values. fig1_construction.py and fig1_construction_v7.py are not changed."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
src = (HERE / "fig1_construction.py").read_text(encoding="utf-8")
EDITS = [('"Same three protein-disjoint partitions in every cell;\\n"\n'
          '                               "each trained ensemble is scored on both test constructions."',
          '"The same three protein-disjoint train–validation–test\\n"\n'
          '                               "splits in every cell; each trained ensemble\\n"\n'
          '                               "is scored on both test constructions."')]
for old, new in EDITS:
    assert src.count(old) == 1, old
    src = src.replace(old, new)
base = {"__name__": "fig1_base", "__file__": str(HERE / "fig1_construction.py")}
exec(compile(src, "fig1_construction", "exec"), base)
INK, INK2, MM, label, save, setup = (base[k] for k in ("INK", "INK2", "MM", "label", "save", "setup"))
T_COL, U_COL = INK2, "#2a78d6"          # as in panel b: threshold reconstruction grey, unrestricted blue


def mini(ax, title, colour, bottom):
    ax.set_xscale("log")
    ax.set_xlim(1, 200)
    ax.set_ylim(0, 1.12)
    ax.set_xticks([1, 10, 100], ["1", "10", "100"])
    ax.set_yticks([0, 1], ["0", "1"])
    ax.set_xlabel("Annotation depth", labelpad=1)
    ax.set_ylabel("Positive fraction\nof retained sites", labelpad=1)
    ax.set_title(title, fontsize=6.8, color=colour, fontweight="bold", pad=4)
    ax.text(0.5, -0.50, bottom, transform=ax.transAxes, ha="center", va="top", fontsize=6.8, color=colour,
            fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def concept(ax_t, ax_u, ax_mid):
    lw = 2.0
    mini(ax_t, "Threshold-sampled training", T_COL, "shallow → positive")
    ax_t.axvspan(1, 10, color=T_COL, alpha=0.08, lw=0)
    ax_t.plot([1, 10], [1, 1], color=T_COL, lw=lw, solid_capstyle="butt")
    ax_t.plot([10, 10], [1, 0.2], color=T_COL, lw=0.9, ls=(0, (2, 1.5)))
    x = np.geomspace(10, 200, 30)
    ax_t.plot(x, 0.2 + 0.15 * np.log10(x / 10) / np.log10(20), color=T_COL, lw=lw)
    ax_t.text(3.2, 0.9, "below the threshold,\nevery retained site\nis positive", ha="center", va="top", fontsize=5.8,
              color=INK2)
    mini(ax_u, "Unrestricted", U_COL, "deeper → more often annotated")
    ax_u.set_ylabel("")
    x = np.geomspace(1, 200, 60)
    ax_u.plot(x, 0.05 + 0.40 * np.log10(x) / np.log10(200), color=U_COL, lw=lw)
    ax_mid.axis("off")
    ax_mid.text(0.5, 0.55, "≠", ha="center", va="center", fontsize=16, color=INK)


def main():
    setup()
    cov = base["coverage"]()
    fig = plt.figure(figsize=(180 * MM, 186 * MM))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.15, 1], height_ratios=[1, 1.05, 0.55], hspace=0.78, wspace=0.42)
    a = fig.add_subplot(gs[0, 0]); base["schematic"](a); label(a, "a", x=-0.02, y=1.0)
    b = fig.add_subplot(gs[0, 1]); base["purepos"](b); label(b, "b", x=-0.30)
    c = fig.add_subplot(gs[1, 0]); base["coverage_panel"](c, cov); label(c, "c", x=-0.22)
    d = fig.add_subplot(gs[1, 1]); base["design"](d); label(d, "d", x=-0.02, y=1.0)
    sub = gs[2, :].subgridspec(1, 5, width_ratios=[0.5, 1, 0.34, 1, 0.6], wspace=0.0)
    et = fig.add_subplot(sub[0, 1]); eu = fig.add_subplot(sub[0, 3]); em = fig.add_subplot(sub[0, 2])
    concept(et, eu, em)
    label(et, "e", x=-0.62, y=1.12)
    save(fig, "Figure1_construction_coverage_v8")


if __name__ == "__main__":
    main()
