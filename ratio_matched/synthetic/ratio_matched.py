#!/usr/bin/env python3
"""Synthetic ratio-matched control (design fixed in PREREG.md).

Three training constructions on the same population, entity split and tests:
  T          negatives only from entities with >= 10 positives (published rule, pool-capped)
  U          negatives from all entities, five per positive (published rule)
  U_matched  the same positives plus exactly T's number of negatives, drawn uniformly
             from every unlabelled item of the training entities

Each is fitted as base / real / perm, once, and scored on the unrestricted and threshold tests.

    python ratio_matched.py            # 45 fits, writes runs.tsv, preds/, summary.txt
"""
import copy
import hashlib
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "analysis"))
from synthetic_reversal import construct, fit_predict, generate, homogeneity, measured_r2  # noqa: E402

THRESHOLD, R2_TARGET, N_ENT, SEEDS = 10, 0.60, 6000, range(5)
TEST_OFFSET, MATCH_OFFSET = 99, 1000
OUT = Path(__file__).resolve().parent
CODE = {"sha_this": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16],
        "sha_synthetic_reversal": hashlib.sha256(
            (REPO / "analysis" / "synthetic_reversal.py").read_bytes()).hexdigest()[:16]}


def population(seed):
    """Population, split and training constructions, in the rng order of one_configuration()."""
    rng = np.random.default_rng(seed)
    d = generate(n_entities=N_ENT, r2_target=R2_TARGET, rng=rng)
    perm = rng.permutation(d["n_entities"])
    test_ent = np.zeros(d["n_entities"], bool)
    test_ent[perm[:d["n_entities"] // 5]] = True
    in_test = test_ent[d["ent"]]

    trains = {}
    for name, t in (("T", THRESHOLD), ("U", 0)):
        sel = construct(d, t, rng=copy.deepcopy(rng))      # both start from the same rng state
        trains[name] = sel[~in_test[sel]]

    # U_matched: T's positives and T's negative count, negatives uniform over training entities
    y = d["y"]
    pos = trains["T"][y[trains["T"]] == 1]
    assert np.array_equal(np.sort(pos), np.sort(trains["U"][y[trains["U"]] == 1]))
    assert np.array_equal(np.sort(pos), np.where((y == 1) & ~in_test)[0])
    n_neg = int((y[trains["T"]] == 0).sum())
    rng_m = np.random.default_rng(seed + MATCH_OFFSET)
    neg = rng_m.choice(np.where((y == 0) & ~in_test)[0], n_neg, replace=False)
    idx = np.concatenate([pos, neg])
    rng_m.shuffle(idx)
    trains["U_matched"] = idx

    tests = {}
    for name, t in (("threshold", THRESHOLD), ("natural", 0)):
        sel = construct(d, t, rng=np.random.default_rng(seed + TEST_OFFSET))
        tests[name] = sel[in_test[sel]]
    return d, trains, tests


def describe(d, tr):
    y, ent = d["y"][tr], d["ent"][tr]
    g = pd.DataFrame({"e": ent, "y": y}).groupby("e")["y"].mean()
    return {"n_pos": int(y.sum()), "n_neg": int((y == 0).sum()), "pos_rate": y.mean(),
            "n_train_entities": len(g), "pure_pos_entity_share": (g == 1.0).mean(),
            "homog_train": homogeneity(d, tr)}


def main():
    rows, t0 = [], time.time()
    (OUT / "preds").mkdir(exist_ok=True)
    for seed in SEEDS:
        d, trains, tests = population(seed)
        r2 = measured_r2(d["feat"], d["depth"], seed)
        te_all = np.concatenate([tests["threshold"], tests["natural"]])
        cut = len(tests["threshold"])
        for cons, tr in trains.items():
            out = {"construction": cons, "seed": seed, "threshold": THRESHOLD, "r2_target": R2_TARGET,
                   "r2_measured": r2, "test_seed": seed + TEST_OFFSET, "init_seed": seed, "perm_seed": seed,
                   "match_seed": seed + MATCH_OFFSET if cons == "U_matched" else np.nan, **CODE,
                   **describe(d, tr)}
            preds = {"te_threshold": tests["threshold"], "te_natural": tests["natural"]}
            for mode in ("base", "real", "perm"):
                p = fit_predict(d, tr, te_all, mode, seed=seed)
                for name, s in (("threshold", p[:cut]), ("natural", p[cut:])):
                    preds[f"{mode}_{name}"] = s
                    yt = d["y"][tests[name]]
                    out[f"auroc_{mode}_{name}"] = roc_auc_score(yt, s)
                    out[f"auprc_{mode}_{name}"] = average_precision_score(yt, s)
                print(f"  seed {seed}  {cons:<9}  {mode:<4}  {time.time() - t0:6.0f} s", flush=True)
            for name in tests:
                for m in ("auroc", "auprc"):
                    b = out[f"{m}_base_{name}"]
                    out[f"d{m}_real_{name}"] = out[f"{m}_real_{name}"] - b
                    out[f"d{m}_perm_{name}"] = out[f"{m}_perm_{name}"] - b
                out[f"real_beats_perm_{name}"] = int(out[f"auroc_real_{name}"] > out[f"auroc_perm_{name}"])
            df = pd.DataFrame({"e": d["ent"][tests["natural"]]})
            for mode in ("base", "real"):
                df["s"] = preds[f"{mode}_natural"]
                g = df.groupby("e")["s"].mean()
                out[f"rho_score_depth_{mode}"] = spearmanr(d["depth"][g.index], g.values).statistic
            out["rho_delta"] = out["rho_score_depth_real"] - out["rho_score_depth_base"]
            np.savez_compressed(OUT / "preds" / f"seed{seed}_{cons}.npz", train_idx=tr, **preds)
            rows.append(out)
    runs = pd.DataFrame(rows)
    runs.to_csv(OUT / "runs.tsv", sep="\t", index=False)
    print(f"done in {time.time() - t0:.0f} s; {len(runs) * 3} fits")


if __name__ == "__main__":
    main()
