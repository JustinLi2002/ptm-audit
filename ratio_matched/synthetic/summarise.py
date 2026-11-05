#!/usr/bin/env python3
"""Reproduction check against the published scan, paired summaries and the decision rule of PREREG.md."""
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
SWEEP = Path(__file__).resolve().parents[2] / "figures" / "figure6_sweep.tsv"
CONS = ["T", "U_matched", "U"]

r = pd.read_csv(OUT / "runs.tsv", sep="\t")
lines = []
say = lines.append

# 1. seeds 0-2 of T and U must reproduce the published scan exactly (same rng sequence)
s = pd.read_csv(SWEEP, sep="\t")
s = s[np.isclose(s.r2_target, 0.60)]
cols = [f"auroc_{m}_natural" for m in ("base", "real", "perm")] + [f"auprc_{m}_natural" for m in ("base", "real", "perm")]
worst = 0.0
for cons, t in (("T", 10), ("U", 0)):
    for seed in (0, 1, 2):
        a = r[(r.construction == cons) & (r.seed == seed)][cols].iloc[0].values
        b = s[(s.threshold == t) & (s.seed == seed)][cols].iloc[0].values
        worst = max(worst, np.abs(a - b).max())
    a = r[(r.construction == "T") & (r.seed.isin([0, 1, 2]))]["auroc_real_threshold"].values
    b = s[(s.threshold == 10)].sort_values("seed")["auroc_real_threshold"].values
    worst = max(worst, np.abs(a - b).max())
say(f"reproduction of figure6_sweep.tsv (seeds 0-2, T and U, natural test; T also threshold test): "
    f"max |diff| = {worst:.2e}")

# 2. per-population table
say("\nper population, unrestricted test (natural):")
key = ["dauroc_real_natural", "dauprc_real_natural", "dauroc_perm_natural", "real_beats_perm_natural",
       "dauroc_real_threshold", "rho_delta"]
desc = ["n_pos", "n_neg", "pos_rate", "n_train_entities", "pure_pos_entity_share", "homog_train"]
for seed in sorted(r.seed.unique()):
    for c in CONS:
        x = r[(r.seed == seed) & (r.construction == c)].iloc[0]
        say(f"  seed {seed} {c:<9} " + " ".join(f"{k.replace('_natural', '')}={x[k]:+.4f}" for k in key)
            + f" | pos_rate={x.pos_rate:.3f} n_neg={x.n_neg} homog={x.homog_train:.2f}%")

say("\nmeans over five populations (min, max):")
g = r.groupby("construction")
for c in CONS:
    x = g.get_group(c)
    say(f"  {c:<9} " + "  ".join(f"{k}={x[k].mean():+.4f} [{x[k].min():+.4f}, {x[k].max():+.4f}]" for k in key[:3])
        + f"  real>perm {int(x.real_beats_perm_natural.sum())}/5")
    say(f"  {'':<9} " + "  ".join(f"{k}={x[k].mean():.3f}" for k in desc) +
        f"  r2_measured={x.r2_measured.mean():.3f}")
    say(f"  {'':<9} absolute AUROC natural: base {x.auroc_base_natural.mean():.4f}, real {x.auroc_real_natural.mean():.4f}, "
        f"perm {x.auroc_perm_natural.mean():.4f}")

# 3. paired differences across populations
p = r.pivot(index="seed", columns="construction", values="dauroc_real_natural")
say("\npaired differences in dAUROC (unrestricted test), per population:")
for a, b in (("U_matched", "T"), ("U", "U_matched"), ("U", "T")):
    dif = p[a] - p[b]
    say(f"  {a} - {b}: " + ", ".join(f"{v:+.4f}" for v in dif) + f"   mean {dif.mean():+.4f}, "
        f"{int((dif > 0).sum())}/5 positive")

# 4. decision rule, in the order fixed in PREREG.md
T, U, M = p["T"], p["U"], p["U_matched"]
if not ((T < 0).all() and (U > 0).all()):
    verdict = "4 (inconclusive): reproduction failed; check implementation"
elif (M > 0).all() and abs(M.mean() - U.mean()) < abs(M.mean() - T.mean()):
    verdict = "2: class ratio and sample size alone do not explain the contrast"
    if (np.sign(M - U) == np.sign((M - U).iloc[0])).all():
        verdict += "; U_matched - U has one sign in all five populations, so matching also changed the effect"
elif M.mean() <= 0 or abs(M.mean() - T.mean()) <= abs(M.mean() - U.mean()):
    verdict = "3: the mechanism attribution must be narrowed"
else:
    verdict = "4 (inconclusive): U_matched of mixed sign"
say(f"\ndecision rule: {verdict}")

txt = "\n".join(lines)
print(txt)
(OUT / "summary.txt").write_text(txt + "\n", encoding="utf-8")
