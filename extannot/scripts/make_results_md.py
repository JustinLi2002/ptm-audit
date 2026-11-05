"""Markdown tables for RESULTS.md, generated from the outputs of extannot_check.py so that no number is copied by hand.

Usage: python -E -P make_results_md.py OUT_DIR > tables.md"""
import sys
from pathlib import Path

import pandas as pd

OUT = Path(sys.argv[1])
TASKS = ["phosphorylation_st", "phosphorylation_y", "acetylation_k", "methylation_k", "methylation_r",
         "sumoylation_k", "ubiquitination_k", "glycosylation_n"]
NAME = {"phosphorylation_st": "Phosphorylation S/T", "phosphorylation_y": "Phosphorylation Y",
        "acetylation_k": "Acetylation K", "methylation_k": "Methylation K/R", "methylation_r": "Methylation R",
        "sumoylation_k": "Sumoylation K", "ubiquitination_k": "Ubiquitination K", "glycosylation_n": "N-glycosylation N"}


def md(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    lines += ["| " + " | ".join(str(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join(lines)


pct = lambda x: f"{100 * x:.2f}"
sg = lambda x: f"{x:+.4f}"

# Step 1: share of negatives (and, for context, positives) annotated in the union; mean over the three splits
c = pd.read_csv(OUT / "step1_counts.tsv", sep="\t").query("source == 'union'")
rows = []
for t in TASKS:
    r = {"Task": NAME[t]}
    for cl in ("T", "U"):
        for pt in ("train", "test"):
            g = c.query("task == @t and construction == @cl and partition == @pt")
            r[f"{cl} {pt} neg %"] = pct(g.neg_share.mean())
    r["test pos %"] = pct(c.query("task == @t and construction == 'U' and partition == 'test'").pos_share.mean())
    d = pd.read_csv(OUT / "step1_test_negatives_by_donor_status.tsv", sep="\t").query("source == 'union' and task == @t")
    for lab in ("donor", "non-donor"):
        r[f"U test neg % ({lab} proteins)"] = pct(d.query("protein_status == @lab").neg_share.mean())
    rows.append(r)
print("### Step 1: share of negatives that are externally annotated (union; mean over three splits)\n")
print(md(pd.DataFrame(rows)) + "\n")

# Step 2
s = pd.read_csv(OUT / "step2_summary.tsv", sep="\t")
rows = []
for fn in ("ESM-2", "ProtT5", "interaction"):
    for src, corr in (("none", "original"), ("union", "relabel"), ("union", "drop")):
        g = s.query("cell == 'TU' and feature == @fn and source == @src and correction == @corr").iloc[0]
        rows.append({"Feature": fn, "Correction": corr, **{NAME[t]: sg(g[t]) for t in TASKS},
                     "8-task mean": sg(g.mean8), "7-unit mean": sg(g.mean7), "Retained": f"{g.retained:.3f}",
                     "Tasks < 0": f"{g.tasks_negative}/8"})
print("### Step 2: TU ΔAUROC (feature − baseline), mean over three splits, union of sources\n")
print(md(pd.DataFrame(rows)) + "\n")

rows = []
for cell in ("TU", "TT", "UU"):
    for src in ("none", "union", "UniProt:primary", "UniProt:primary+similarity", "CPLM:primary",
                "PhosphoSitePlus:primary", "PhosphoSitePlus:low_throughput"):
        if src not in set(s.source):
            continue
        for corr in (("original",) if src == "none" else ("relabel", "drop")):
            g = s.query("cell == @cell and feature == 'ESM-2' and source == @src and correction == @corr").iloc[0]
            rows.append({"Cell": cell, "Source": "—" if src == "none" else src, "Correction": corr,
                         "Negatives flagged (3 splits)": g.n_flagged, "8-task mean": sg(g.mean8),
                         "7-unit mean": sg(g.mean7), "Retained": f"{g.retained:.3f}",
                         "Tasks < 0": f"{g.tasks_negative}/8"})
print("### Step 2: ESM-2 by cell and source\n")
print(md(pd.DataFrame(rows)) + "\n")

# Enrichment
e = pd.read_csv(OUT / "enrichment_summary.tsv", sep="\t").query("source == 'union'").set_index("task")
rows = [{"Task": NAME[t], "Baseline": f"{e.loc[t, 'enrich_baseline']:.3f}", "ESM-2": f"{e.loc[t, 'enrich_ESM-2']:.3f}",
         "ESM-2 − baseline": sg(e.loc[t, "ESM-2_minus_baseline"]), "ProtT5": f"{e.loc[t, 'enrich_ProtT5']:.3f}",
         "Interaction": f"{e.loc[t, 'enrich_interaction']:.3f}"} for t in TASKS]
print("### Enrichment: AUROC of the TU score for annotated vs other unrestricted test negatives (union; mean over "
      "three splits)\n")
print(md(pd.DataFrame(rows)))
