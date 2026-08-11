#!/usr/bin/env python3
"""Patch make_figure_merged.py: feature-family detection and legend encoding.

Two changes.

1. load() classified any file not ending in '__esm' as the interaction
   embedding.  That was correct while only two families existed; with the
   ProtT5 runs in the same directory it silently overwrites the interaction
   series with ProtT5 values.  Files are now matched against an explicit list
   of suffixed families, longest suffix first.

2. The legend used filled circles for the two feature families, the same glyph
   the plot uses for "real feature beats its permutation", so the colour and
   shape channels were not separable by eye.  Worse, the handle for "beats its
   permutation" was drawn hollow while the plot draws it filled, so the legend
   contradicted the figure.  Families now use colour swatches and the
   win/lose handles are neutral grey, which reads as "shape only; colour is
   carried separately".
"""
from pathlib import Path

P = Path("make_figure_merged.py")
s = P.read_text(encoding="utf-8")
log = []


def sub(old, new, label):
    global s
    n = s.count(old)
    assert n == 1, f"[{label}] 期望 1 处,实际 {n}"
    s = s.replace(old, new, 1)
    log.append(label)


# ── 1. 特征族识别 ───────────────────────────────────────────────────
sub("from matplotlib.lines import Line2D",
    "from matplotlib.lines import Line2D\nfrom matplotlib.patches import Patch",
    "import Patch")

sub("C_PPI, C_ESM, C_PART = '#1f77b4', '#d62728', '#bbbbbb'",
    "C_PPI, C_ESM, C_PART = '#1f77b4', '#d62728', '#bbbbbb'\n"
    "C_MARK = '0.30'          # shape channel only; colour means the family\n"
    "# Suffixed feature families. 'ppi' carries no suffix and holds the\n"
    "# baselines, so anything unmatched must fall through to it -- which is\n"
    "# why the list has to be explicit rather than a single '__esm' test.\n"
    "SUFFIXED = ('prott5', 'esm')",
    "定义 SUFFIXED 与 C_MARK")

sub("""        feat = 'esm' if os.path.basename(f)[:-5].endswith('__esm') else 'ppi'""",
    """        stem = os.path.basename(f)[:-5]
        feat = next((k for k in SUFFIXED if stem.endswith('__' + k)), 'ppi')""",
    "load() 按后缀表识别特征族")

# ── 2. 图例 ─────────────────────────────────────────────────────────
sub("""    fig.legend(handles=[
        Line2D([], [], marker='o', ls='', color=C_PPI, label='interaction embedding'),
        Line2D([], [], marker='o', ls='', color=C_ESM, label='frozen ESM-2'),
        Line2D([], [], marker='o', ls='', mfc='none', mec='k', color='k',
               label='real feature beats its permutation'),
        Line2D([], [], marker='v', ls='', mfc='none', mec='k', color='k',
               label='permutation beats the real feature'),
        Line2D([], [], marker='o', ls='', color=C_PART, label='individual partition'),
    ], frameon=False, fontsize=6.5, loc='lower center', ncol=3,
        bbox_to_anchor=(0.5, -0.06))""",
    """    fig.legend(handles=[
        Patch(facecolor=C_PPI, edgecolor='none', label='interaction embedding'),
        Patch(facecolor=C_ESM, edgecolor='none', label='frozen ESM-2'),
        Line2D([], [], marker='o', ls='', mfc=C_MARK, mec=C_MARK,
               label='real feature beats its permutation'),
        Line2D([], [], marker='v', ls='', mfc='none', mec=C_MARK,
               label='permutation beats the real feature'),
        Line2D([], [], marker='.', ls='', color=C_PART,
               label='individual partition'),
    ], frameon=False, fontsize=6.5, loc='lower center', ncol=3,
        bbox_to_anchor=(0.5, -0.06))""",
    "图例:特征族用色块,胜负用中性灰形状")

P.write_text(s, encoding="utf-8")
print(f"完成 {len(log)} 处:")
for i, t in enumerate(log, 1):
    print(f"  {i}. {t}")
