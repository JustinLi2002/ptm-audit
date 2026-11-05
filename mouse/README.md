# Independent mouse dataset

The eight human tasks come from one benchmark and one organism, so they are not independent tests of the effect of the
donor rule. This folder applies the same negative-sampling intervention to phosphorylation S/T sites from an independent
experimental source: Supplementary Table 2 of Giansanti et al., Nat Methods 2022;19:803-11,
doi:10.1038/s41592-022-01526-y, whose data are deposited in PRIDE under accession PXD030983.

The ten-site threshold is our intervention on these data. We do not suggest that the original study used such a rule.

## Result

ESM-2 only. Decision row B of the rule fixed in `PREREG_mouse.md`; matched control row 2. All six cells are in
`results/mouse_results.tsv` and Supplementary Table 17.

This is a replication of the training-construction effect on an independent source. It is not a full replication of the
human pattern: a benefit on the threshold-sampled test was not established, because the TT interval includes zero, and
the data do not show that the benefit is absent either.

`LABEL_CORRECTION_rowB_20261103.md` records a correction to the wording of the row B label, made after the results were
known. The condition, the observed row and every analysis are unchanged, and the text of `PREREG_mouse.md` is unchanged.

## Why the constructions are archived

`rebuild_datasets.py` enumerates candidate negatives by iterating over a Python set of accessions, so the order depends
on string hashing. Rerunning the scripts on the source table does not reproduce the sampled negatives. The files actually
used are therefore in the data record, not only the scripts and the source.

## Contents

| Path | What it is |
|---|---|
| `PREREG_mouse.md`, `PREREG_mouse.sha256` | Design, outcomes and decision rule, fixed before any site was mapped |
| `RUNLOG.md` | Every step, in order |
| `LABEL_CORRECTION_rowB_20261103.md` | The wording correction described above |
| `scripts/map_sites.py` | Source table to Swiss-Prot sites; deterministic, no sampling |
| `scripts/analyse_mouse.py` | Outcomes and decision rows |
| `scripts/hrp/` | Cluster jobs: matched sets, training (`train_mouse.py`), ESM-2 extraction, splits |
| `scripts/make_bundle.py`, `submit_template.sh`, `fetch_*.sh` | How the run was submitted and fetched |
| `results/` | Tables, construction prediction, mapping report, checksums, scheduler accounting |

## Reproducing the numbers

Extract `mouse_predictions.tar.gz` and `mouse_constructions.tar.gz` from the data record (concept DOI
10.5281/zenodo.21670042) in one directory, then run from the repository root:

    python mouse/scripts/analyse_mouse.py <dir>/mouse <output dir>

No training is needed. `build_matched_mouse.py`, `train_mouse.py` and the slurm files need `PTM_AUDIT_BASE` set to a
directory named `mouse_v1`; this guard kept the run from writing to the human tree.

## Source table and licence

The source workbook is not redistributed; obtain it from the publisher of the article above. Its sha256 is in
`results/mapping_report.json`. The article carries no open licence that we could find. The data record holds only the
protein, position and residue of each site after mapping to UniProtKB/Swiss-Prot, and how it mapped (`sites_mapped.tsv`,
and the labels in the construction and prediction files), so that the analysis can be checked. See the record's README for attribution and
the removal clause.
