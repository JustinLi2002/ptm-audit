#!/usr/bin/env python3
"""Real-data ratio-matched control: per-partition contributions and the decision rule of PREREG_real.md.

Reads the fetched run summaries (fetch/, HRP-relative paths) and writes runs_real.tsv and summary_real.txt.

    python analyse_real.py [DIR]

DIR holds pdisjoint_runs_v2 and pdisjoint_runs_matched_v4 from the data record; the default is fetch/ beside this script.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FETCH = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "fetch"
CROSS = Path(__file__).resolve().parents[2] / "results" / "crosseval_verify.tsv"
TASKS = ["acetylation_k", "phosphorylation_st"]
CONS = {"T": ("pdisjoint_runs_v2", "replica"), "matched": ("pdisjoint_runs_matched_v4", "matched"),
        "U": ("pdisjoint_runs_v2", "rebuilt")}
TESTS = {"unrestricted": "rebuilt", "threshold": "replica"}


def run(folder, ds, task, cond, s):
    suffix = "" if cond == "baseline" else "__esm"
    p = FETCH / folder / f"{task}__{ds}__{cond}__split{s}{suffix}.json"
    return json.loads(p.read_text()) if p.exists() else None


rows, missing = [], []
for task in TASKS:
    for cons, (folder, ds) in CONS.items():
        for s in (0, 1, 2):
            r = {c: run(folder, ds, task, c, s) for c in ("baseline", "ppi", "shuffled")}
            if any(v is None for v in r.values()):
                missing.append((task, cons, s, [c for c, v in r.items() if v is None]))
                continue
            out = dict(task=task, construction=cons, split=s, pos_train=r["baseline"]["pos_train"],
                       n_train=r["baseline"]["n_train"])
            for tn, key in TESTS.items():
                b, e, q = (r[c]["tests"][key] for c in ("baseline", "ppi", "shuffled"))
                out[f"auroc_base_{tn}"] = b["auroc"]
                out[f"dauroc_{tn}"] = e["auroc"] - b["auroc"]
                out[f"dauprc_{tn}"] = e["auprc"] - b["auprc"]
                out[f"dauroc_perm_{tn}"] = q["auroc"] - b["auroc"]
                out[f"real_beats_perm_{tn}"] = int(e["auroc"] > q["auroc"])
            rows.append(out)

lines = []
say = lines.append
if missing:
    say("MISSING RUNS: " + "; ".join(f"{t} {c} split{s} {m}" for t, c, s, m in missing))
d = pd.DataFrame(rows)
d.to_csv(HERE / "runs_real.tsv", sep="\t", index=False)

# the existing T and U values must agree with the manuscript's cross-evaluation table
v = pd.read_csv(CROSS, sep="\t").drop_duplicates(["task", "train", "eval", "split", "family"])
v = v[(v.family == "esm2") & (v["eval"] == "rebuilt")].set_index(["task", "train", "split"]).d_auroc
dev = [abs(x.dauroc_unrestricted - v[(x.task, {"T": "replica", "U": "rebuilt"}[x.construction], x.split)])
       for x in d[d.construction != "matched"].itertuples()]
say(f"T/U dAUROC vs crosseval_verify.tsv: max |diff| = {max(dev):.4f} over {len(dev)} task-partitions")

for task in TASKS:
    say(f"\n== {task}")
    t = d[d.task == task]
    for cons in CONS:
        x = t[t.construction == cons].sort_values("split")
        if x.empty:
            continue
        say(f"  {cons:<8} pos_train {x.pos_train.mean():.3f}  n_train {x.n_train.mean():,.0f}")
        for col in ("dauroc_unrestricted", "dauprc_unrestricted", "dauroc_perm_unrestricted",
                    "dauroc_threshold", "dauprc_threshold"):
            say(f"           {col:<26} mean {x[col].mean():+.4f}  per partition "
                + ", ".join(f"{val:+.4f}" for val in x[col]))
        say(f"           real > permuted (unrestricted test): {int(x.real_beats_perm_unrestricted.sum())}/{len(x)}"
            f"; baseline AUROC {x.auroc_base_unrestricted.mean():.4f}")
    p = t.pivot(index="split", columns="construction", values="dauroc_unrestricted")
    if not {"T", "U", "matched"} <= set(p.columns) or p.isna().any().any():
        say("  decision rule: not applied (runs missing)")
        continue
    T, U, M = p["T"], p["U"], p["matched"]
    say("  paired, per partition: matched - T " + ", ".join(f"{z:+.4f}" for z in M - T)
        + " | U - matched " + ", ".join(f"{z:+.4f}" for z in U - M))
    if not (T.mean() < 0 and U.mean() > 0):
        verdict = "1: this task cannot address the question (T not negative or U not positive)"
    elif (M > 0).all() and abs(M.mean() - U.mean()) < abs(M.mean() - T.mean()):
        verdict = "2: class ratio and training-set size alone do not explain the contrast"
    elif M.mean() <= 0 or abs(M.mean() - T.mean()) <= abs(M.mean() - U.mean()):
        verdict = "3: the mechanism attribution must be narrowed for this task"
    else:
        verdict = "4: mixed; report effect sizes"
    say(f"  decision rule: {verdict}")

txt = "\n".join(lines)
print(txt)
(HERE / "summary_real.txt").write_text(txt + "\n", encoding="utf-8")
