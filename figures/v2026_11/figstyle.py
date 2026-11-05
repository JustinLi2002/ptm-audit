"""Shared style and data paths for the redesigned figures (draft for the restructure)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
RES = REPO / "results"
S3 = RES / "section3_20260920"
OUT = Path(__file__).resolve().parent

# ink and chrome (light surface only: these are print figures)
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
# channel identity: the first three categorical slots, validated all-pairs on white
CH = {"interaction": "#2a78d6", "esm2": "#eb6834", "prott5": "#1baf7a"}
CH_LABEL = {"interaction": "Interaction", "esm2": "ESM-2", "prott5": "ProtT5"}
# ordinal ramp for donor threshold T (validated --ordinal on white)
RAMP5 = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]

TASKS = ["phosphorylation_st", "phosphorylation_y", "acetylation_k", "methylation_k",
         "methylation_r", "sumoylation_k", "ubiquitination_k", "glycosylation_n"]
SHORT = {"phosphorylation_st": "Phospho S/T", "phosphorylation_y": "Phospho Y",
         "acetylation_k": "Acetyl K", "methylation_k": "Methyl K/R",
         "methylation_r": "Methyl R", "sumoylation_k": "Sumo K",
         "ubiquitination_k": "Ubiq K", "glycosylation_n": "N-Glyc N"}
# the four cells of the design: (training construction, evaluation construction)
CELLS = [("replica", "replica"), ("replica", "rebuilt"), ("rebuilt", "replica"), ("rebuilt", "rebuilt")]
CELL_NAME = {("replica", "replica"): "Threshold train · threshold test",
             ("replica", "rebuilt"): "Threshold train · unrestricted test",
             ("rebuilt", "replica"): "Unrestricted train · threshold test",
             ("rebuilt", "rebuilt"): "Unrestricted train · unrestricted test"}

MM = 1 / 25.4


def setup():
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 7, "axes.titlesize": 7.5,
        "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5, "axes.edgecolor": AXIS, "axes.linewidth": 0.6,
        "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.6,
        "legend.frameon": False, "savefig.facecolor": "white",
        "figure.facecolor": "white", "lines.linewidth": 1.4,
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })


def label(ax, letter, x=-0.12, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left", color=INK)


def hgrid(ax):
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def vgrid(ax):
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
