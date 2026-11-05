#!/usr/bin/env python3
"""Map the Giansanti et al. (2022) organ phosphosites to UniProtKB/Swiss-Prot mouse canonical sequences
(PREREG_mouse.md, "Data").

Site filters (as in the preflight): localization probability >= 0.75, an intensity in at least one of the 41
tissues, a single leading protein, and the centre of the 31-residue window equal to the reported residue.

Each site is placed by its window, not by its 2018 position. The run of residues in the window (with "_" padding
read as a sequence end) must occur exactly once in the Swiss-Prot canonical entry of an accession in the row's
"Proteins" field (isoform accessions count through their canonical entry). The leading accession's entry is tried
first; otherwise exactly one of the other entries must match. Anything else is dropped and counted.

Writes, in OUT_DIR:
  sites_mapped.tsv                     one row per mapped site (accession, pos, aa, how)
  retrain/phosphorylation_{st,y}_{train,test}.tsv   positives in the layout rebuild_datasets.py reads
  reference.fasta                      Swiss-Prot entries carrying a mapped site, UniProt headers kept
  mapping_report.json

Usage: python map_sites.py WORKBOOK SWISSPROT_FASTA OUT_DIR"""
import hashlib
import json
import re
import sys
from pathlib import Path

import pandas as pd

WB, FA, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
SHEET = "Supplementary Table 2 - Organs"
HALF = 15


def read_fasta(path):
    seqs, hdrs, acc, buf = {}, {}, None, []
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if acc:
                seqs[acc] = "".join(buf)
            acc, buf = line.split("|")[1], []
            hdrs[acc] = line
        else:
            buf.append(line.strip())
    seqs[acc] = "".join(buf)
    return seqs, hdrs


def locate(seq, win):
    """1-based positions in seq at which the window's centre falls; [] if none."""
    left, right = win[:HALF].lstrip("_"), win[HALF + 1:].rstrip("_")
    n_start, c_end = win.startswith("_"), win.endswith("_")
    core = left + win[HALF] + right
    hits, i = [], seq.find(core)
    while i >= 0:
        if (not n_start or i == 0) and (not c_end or i + len(core) == len(seq)):
            hits.append(i + len(left) + 1)
        i = seq.find(core, i + 1)
    return hits


sp, hdr = read_fasta(FA)
df = pd.read_excel(WB, sheet_name=SHEET, engine="openpyxl")
inten = [c for c in df.columns if c.startswith("Intensity (Log10)")]
assert len(inten) == 41, len(inten)
rep = {"rows": len(df)}
d = df[df["Localization prob"] >= 0.75]
d = d[d[inten].notna().any(axis=1)]
rep["localized_with_tissue_intensity"] = len(d)
d = d[~d["Leading proteins"].astype(str).str.contains(";")]
rep["single_leading_protein"] = len(d)
d = d[d["Sequence window"].str[HALF] == d["Amino acid"]]
rep["window_centre_agrees"] = len(d)
d = d.drop_duplicates(["Leading proteins", "Position", "Amino acid"])
rep["unique_sites_before_mapping"] = len(d)
rep["unique_ST_before_mapping"] = int(d["Amino acid"].isin(["S", "T"]).sum())

base = lambda a: re.sub(r"-\d+$", "", a)
rows, why = [], {"leading": 0, "other_entry": 0, "no_swissprot_entry": 0, "no_match": 0, "ambiguous": 0}
for r in d.itertuples(index=False):
    win, aa = r[6], r[7]
    lead = base(str(r[2]))
    cands = [lead] + [a for a in dict.fromkeys(base(x) for x in str(r[0]).split(";")) if a != lead]
    cands = [a for a in cands if a in sp]
    if not cands:
        why["no_swissprot_entry"] += 1
        continue
    hits = {a: locate(sp[a], win) for a in cands}
    if lead in hits and len(hits[lead]) == 1:
        acc, pos, how = lead, hits[lead][0], "leading"
    else:
        ok = [(a, h[0]) for a, h in hits.items() if len(h) == 1 and a != lead]
        many = [a for a, h in hits.items() if len(h) > 1]
        if len(ok) == 1 and not many:
            (acc, pos), how = ok[0], "other_entry"
        else:
            why["ambiguous" if ok or many else "no_match"] += 1
            continue
    assert sp[acc][pos - 1] == aa
    why[how] += 1
    rows.append((acc, pos, aa, how))

m = pd.DataFrame(rows, columns=["accession", "pos", "aa", "how"]).drop_duplicates(["accession", "pos"])
rep["mapping"] = why
rep["mapped_unique_sites"] = len(m)
st = m[m.aa.isin(["S", "T"])]
rep["mapped_ST_sites"] = len(st)
rep["mapped_ST_share_of_prefiltered_ST"] = len(st) / rep["unique_ST_before_mapping"]
depth = m.groupby("accession").size()
donors = set(depth.index[depth >= 10])
rep["ST_proteins"] = int(st.accession.nunique())
rep["ST_proteins_merged_STY_depth_ge10"] = int(st.accession.isin(donors).groupby(st.accession).first().sum())
rep["ST_sites_on_depth_ge10_proteins"] = int(st.accession.isin(donors).sum())
rep["ST_sites_on_depth_lt10_proteins"] = int((~st.accession.isin(donors)).sum())
rep["Y_sites"] = int((m.aa == "Y").sum())
rep["stop_rule_80pct_met"] = rep["mapped_ST_share_of_prefiltered_ST"] >= 0.80

OUT.mkdir(parents=True, exist_ok=True)
m.to_csv(OUT / "sites_mapped.tsv", sep="\t", index=False)
(OUT / "retrain").mkdir(exist_ok=True)
for task, res in (("phosphorylation_st", "ST"), ("phosphorylation_y", "Y")):
    x = m[m.aa.isin(list(res))].rename(columns={"accession": "protein"})[["protein", "aa", "pos"]].assign(y=1)
    x.to_csv(OUT / "retrain" / f"{task}_train.tsv", sep="\t", index=False)
    x.iloc[:0].to_csv(OUT / "retrain" / f"{task}_test.tsv", sep="\t", index=False)
with open(OUT / "reference.fasta", "w", encoding="utf-8", newline="\n") as fh:
    for a in sorted(m.accession.unique()):
        fh.write(hdr[a] + "\n")
        s = sp[a]
        fh.write("\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
rep["inputs"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (WB, FA)}
(OUT / "mapping_report.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")
print(json.dumps(rep, indent=2))
