"""List the unrestricted test negatives flagged by the union of sources (primary tiers), one row per task, split,
protein and position, with the same matching rules as extannot_check.py (accession present in the reference,
residue at the position equal to the database residue, modification kind and residue of the task). Used to report
how many flagged observations are distinct sites, since test partitions of different splits overlap.

Usage: python -X utf8 -E -P list_flagged.py ZENODO_X SITES_TSV OUT_TSV"""
import sys
from pathlib import Path

import pandas as pd

ZX, SITES, OUT = map(Path, sys.argv[1:4])
TASKS = ["phosphorylation_st", "phosphorylation_y", "acetylation_k", "methylation_k", "methylation_r",
         "sumoylation_k", "ubiquitination_k", "glycosylation_n"]
KINDS = {"phosphorylation_st": {"phospho"}, "phosphorylation_y": {"phospho"}, "acetylation_k": {"acetyl"},
         "methylation_k": {"methyl_K", "methyl_R"}, "methylation_r": {"methyl_R"}, "sumoylation_k": {"sumo"},
         "ubiquitination_k": {"ubiquitin"}, "glycosylation_n": {"nglyc"}}
RESIDUES = {"phosphorylation_st": "ST", "phosphorylation_y": "Y", "acetylation_k": "K", "methylation_k": "KR",
            "methylation_r": "R", "sumoylation_k": "K", "ubiquitination_k": "K", "glycosylation_n": "N"}

ref, acc = {}, None
for line in open(ZX / "pdisjoint" / "proteins.fasta", encoding="utf-8"):
    line = line.strip()
    if line.startswith(">"):
        acc = line[1:].split()[0]
        ref[acc] = []
    elif acc:
        ref[acc].append(line)
ref = {k: "".join(v) for k, v in ref.items()}

sites = pd.read_csv(SITES, sep="\t", keep_default_na=False).query("tier == 'primary'")
hit = {t: set() for t in TASKS}
for r in sites.itertuples(index=False):
    seq = ref.get(r.accession)
    if seq is None or not (0 < r.position <= len(seq)) or seq[r.position - 1] != r.residue:
        continue
    for t in TASKS:
        if r.kind in KINDS[t] and r.residue in RESIDUES[t]:
            hit[t].add((r.accession, r.position))

rows = []
for t in TASKS:
    for sp in (0, 1, 2):
        d = pd.read_csv(ZX / "pdisjoint_runs_v2" / f"{t}__replica__baseline__split{sp}__on_rebuilt.pred.tsv",
                        sep="\t", usecols=["protein", "pos", "y"])
        for p, pos in zip(d.protein[d.y == 0], d.pos[d.y == 0]):
            if (p, pos) in hit[t]:
                rows.append((t, sp, p, pos))
out = pd.DataFrame(rows, columns=["task", "split", "protein", "pos"])
out.to_csv(OUT, sep="\t", index=False)
print(f"{len(out)} flagged test-negative observations; "
      f"{len(out[['task', 'protein', 'pos']].drop_duplicates())} distinct task-protein-position combinations")
