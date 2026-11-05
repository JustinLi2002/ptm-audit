"""Figure 3 (v4). As fig3_explanation.py without the across-task panel (old g, now Supplementary Fig. 3),
with the remaining six panels enlarged."""
import matplotlib.pyplot as plt

from fig3_explanation import depth_bins, label_rate, neighbours, probe_r2, probe_scatter, score_shift
from figstyle import MM, label, save, setup


def main():
    setup()
    fig = plt.figure(figsize=(180 * MM, 150 * MM))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.1], hspace=0.5, wspace=0.6)
    a = fig.add_subplot(gs[0, 0]); label_rate(a); label(a, "a", x=-0.3, y=1.06)
    b = fig.add_subplot(gs[0, 1]); probe_scatter(b); label(b, "b", x=-0.3, y=1.06)
    c = fig.add_subplot(gs[0, 2]); probe_r2(c); label(c, "c", x=-0.55, y=1.06)
    d = fig.add_subplot(gs[1, 0]); score_shift(d); label(d, "d", x=-0.42, y=1.06)
    e = fig.add_subplot(gs[1, 1]); neighbours(e); label(e, "e", x=-0.42, y=1.06)
    f = fig.add_subplot(gs[1, 2]); depth_bins(f); label(f, "f", x=-0.3, y=1.06)
    save(fig, "Figure3_explanation_v4")


if __name__ == "__main__":
    main()
