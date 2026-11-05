# External-annotation check: results with UniProt, CPLM and PhosphoSitePlus phosphorylation

Run on 2026-11-02 under PREREG.md. This run supersedes the earlier UniProt + CPLM run in
`results_uniprot_cplm_20261102/`, which is kept for reference. Only one PhosphoSitePlus dataset, phosphorylation,
could be used. The other four could not be obtained (see "Deviations").

## Primary analysis (TU, ESM-2, relabel, union of sources): harm survives

| | Original | Corrected |
|---|---|---|
| 8-task mean ΔAUROC | −0.1008 | −0.0986 |
| 7-unit mean ΔAUROC | −0.0814 | −0.0791 |
| Retained share | — | 0.972 |
| Tasks below zero | 7/8 | 7/8 |

Under the decision rule (mean below zero, retained share ≥ 0.5, at least six of eight tasks below zero), the harm
survives. Step 3, retraining on cleaned negatives, is not triggered. Adding PhosphoSitePlus moved the corrected mean
from −0.0989 (UniProt + CPLM) to −0.0986.

## Inputs

- **Zenodo data record 21670043:** rebuilt_datasets, replica_reconstruction, protein_disjoint_splits and
  crosseval_predictions_v2. All four MD5 checksums match the record.
- **UniProtKB/Swiss-Prot:** release 2026_03, reference proteome UP000005640, reviewed entries. Fetched on 2026-11-02
  through the REST search endpoint: 20,416 entries, with the modified-residue, glycosylation and cross-link features.
- **CPLM:** Homo_sapiens.zip (16,739,171 bytes), downloaded on 2026-11-02. The file inside is dated 2021-06-13.
- **PhosphoSitePlus:** Phosphorylation_site_dataset.
  - The copy J.L. downloaded on 2025-06-06 for an earlier project (`Downloads/01-contextPTM/`).
  - Stored here as `data/psp/Phosphorylation_site_dataset_dl20250606` (25,883,472 bytes, MD5
    585d8cecaa587b80f545d68d136db617).
  - Contents: 334,773 sites, of which 211,139 are human and 11,937 of those have LT_LIT ≥ 1. The file is complete:
    every row has 14 columns, and the last gene is ZZZ3.
- **Reference sequences:** pdisjoint/proteins.fasta from the data record (16,336 proteins).

## Checks made before any correction

- **AUROC reproduction.** The deposited ensemble AUROCs were reproduced from the per-site `y_pred` columns, with a
  maximum absolute difference of 2.2e-16. The uncorrected means reproduce the manuscript: TU −0.1008, TT +0.0437,
  UU +0.0376.
- **Data consistency.**
  - The residue of every benchmark row agrees with the reference sequence.
  - The positives are identical in T and U.
  - Every threshold-sampled negative lies on a donor protein (at least ten positives of the merged modification type),
    which confirms that the donor rule was reimplemented correctly.
- **Matching the external sites.** Accession, position and residue were required to match.

  | Source | Matched | Dropped (residue mismatch) | Outside the reference |
  |---|---|---|---|
  | UniProt, experimental evidence | 42,366 | 497 | 259 |
  | CPLM | 220,617 | 443 | 35,847 |
  | PhosphoSitePlus | 192,562 | 2,475 | 9,144 |

  For PhosphoSitePlus, 6,958 further rows with isoform accessions were not used.
- **Sequence drift between sources (diagnostic, post hoc).** For each CPLM site that matched on residue and has a
  full 15-residue window, the window in CPLM's own sequence was compared with ours. 215,816 windows are identical and
  60 differ, so sequence drift between the sources is negligible.

## What the results show

1. **The correction barely moves the harm.**
   - Across the three splits, 36,828 TU test negatives are externally annotated.
   - Relabelling them as positives moves the ESM-2 eight-task mean from −0.1008 to −0.0986.
   - Dropping them makes the mean slightly more negative (−0.1019).
   - ProtT5 and the interaction embedding behave the same way, with retained shares of 0.971 and 0.963.
   - Only one cell changes sign: the interaction embedding on phosphorylation Y, from −0.0001 to +0.0056 under
     relabelling.
2. **Phosphorylation: PhosphoSitePlus raises the annotated share of the negatives by more than twentyfold, and the
   result does not move.**
   - With UniProt alone, 0.25–0.45% of phosphorylation negatives were annotated. With PhosphoSitePlus, 6.3% of
     unrestricted S/T test negatives and 14.8% of Y test negatives are annotated.
   - Even so, ESM-2 TU ΔAUROC moves only from −0.0266 to −0.0247 in S/T, and from +0.0149 to +0.0160 in Y.
3. **The external annotation is denser in the threshold-sampled negatives in seven of eight tasks.**
   - The share of annotated test negatives is higher in T than in U in every task except phosphorylation Y. In the
     tasks where the shares are large, the figures are, as T test against U test:
     - acetylation K: 15.4% against 8.2%;
     - sumoylation K: 8.5% against 4.0%;
     - ubiquitination K: 15.5% against 10.1%;
     - phosphorylation S/T: 7.3% against 6.3%.
   - Phosphorylation Y is the exception: 13.6% in T against 14.8% in U.
   - Within the unrestricted test negatives, those on donor proteins are annotated at a higher rate than those on
     non-donor proteins in every task:
     - about 1.4 to 1.5 times in the two phosphorylation tasks;
     - 2.3 to 3.2 times in acetylation, sumoylation, ubiquitination and N-glycosylation;
     - more in the two methylation tasks, at low absolute rates.
   - So annotated non-sites are recorded at higher rates on well-studied proteins, which are the proteins the threshold selects. The Methods
     sentence that the residual noise "need not be equal" between constructions can now be given a direction and a
     size.
4. **The direct form of the alternative explanation runs the wrong way.**
   - If the feature pushed true but unannotated sites up, the annotated negatives would rank higher under ESM-2 than
     under the baseline.
   - In six of eight tasks they rank lower. For example, the enrichment AUROC falls from 0.668 to 0.534 in
     acetylation K and from 0.578 to 0.420 in methylation R.
   - ESM-2 ranks them higher only in phosphorylation Y (+0.016), the one task where TU ESM-2 is not harmful, and
     marginally in ubiquitination K (+0.007).
   - This pattern is consistent with the shortcut described in the paper: a threshold-trained feature learns that
     negatives lie on donor proteins, scores negatives on well-studied proteins low, and the annotated negatives
     are recorded at higher rates on exactly those proteins. This reading is an interpretation, not a test.
5. **The in-construction cells are hardly affected.**
   - TT: relabelling lowers the ESM-2 gain from +0.0437 to +0.0389 (retained share 0.866), and all eight tasks stay
     positive. A small part of the gain on threshold-sampled test sets comes from ranking annotated negatives low.
   - UU: the gain is unchanged (+0.0377 after relabelling).

## Deviations from PREREG.md and additions to it

- **PhosphoSitePlus coverage.** Only the phosphorylation dataset was used, in its June 2025 release. The acetylation,
  methylation, ubiquitination and sumoylation datasets could not be obtained:
  - The site relaunched as version 7.0 on 2026-10-01, and the downloads view of the new site is switched off.
  - The current Terms and Conditions state that the licence gives no right to download data from the site without a
    separate written agreement with CST.
  - For the lysine tasks, CPLM is therefore the main source. Methylation R and N-glycosylation rely on UniProt alone.
- **Release date.** The PhosphoSitePlus file is dated by its download (2025-06-06), not at download as PREREG
  anticipated. Its release date is not recorded in the file.
- **Python flags.** Python was run with `-X utf8 -E -P` rather than `-I`, because numpy is installed in the user
  site-packages. With these flags, neither the script directory nor the working directory is on the import path. This
  has no effect on the results.
- **Additions (descriptive, not pre-registered):**
  - validation-partition counts and the annotated share of the positives in step1_counts.tsv;
  - the CPLM window diagnostic;
  - UniProt similarity-only evidence as a separate source set;
  - ProtT5 and the interaction embedding in the enrichment analysis, which PREREG specified for the baseline and
    ESM-2 only.
- **CPLM evidence.** CPLM evidence includes "Dat." (integrated from other databases). As PREREG specifies, all sites
  were used. No experimental-only CPLM tier was run; one would be labelled post hoc.

## Files in this folder

- `report.txt`: checks, matching counts and the primary verdict.
- `sites_report.txt`, `sites.tsv`: the parsed external sites.
- `step1_counts.tsv`: counts by task, split, construction, partition and source.
- `step1_test_negatives_by_donor_status.tsv`: unrestricted test negatives on donor and non-donor proteins.
- `step2_rescoring_by_split.tsv`, `step2_summary.tsv`: the rescoring, per split and summarised.
- `enrichment_by_split.tsv`, `enrichment_summary.tsv`: the enrichment analysis, per split and summarised.

The scripts are in `../scripts/`.

## Tables (generated from the TSV files by make_results_md.py)

### Step 1: share of negatives that are externally annotated (union; mean over three splits)

| Task | T train neg % | T test neg % | U train neg % | U test neg % | test pos % | U test neg % (donor proteins) | U test neg % (non-donor proteins) |
|---|---|---|---|---|---|---|---|
| Phosphorylation S/T | 7.30 | 7.26 | 6.33 | 6.33 | 48.36 | 7.24 | 5.19 |
| Phosphorylation Y | 13.56 | 13.56 | 14.84 | 14.84 | 70.69 | 16.06 | 10.96 |
| Acetylation K | 14.49 | 15.42 | 8.11 | 8.19 | 73.24 | 15.16 | 6.61 |
| Methylation K/R | 1.46 | 0.85 | 0.46 | 0.37 | 13.61 | 1.33 | 0.31 |
| Methylation R | 0.25 | 0.34 | 0.07 | 0.09 | 5.87 | 0.54 | 0.05 |
| Sumoylation K | 8.38 | 8.53 | 3.80 | 3.96 | 71.98 | 8.31 | 2.58 |
| Ubiquitination K | 15.59 | 15.51 | 10.14 | 10.05 | 79.67 | 15.51 | 6.64 |
| N-glycosylation N | 0.73 | 0.91 | 0.68 | 0.45 | 25.31 | 1.01 | 0.39 |

### Step 2: TU ΔAUROC (feature − baseline), mean over three splits, union of sources

| Feature | Correction | Phosphorylation S/T | Phosphorylation Y | Acetylation K | Methylation K/R | Methylation R | Sumoylation K | Ubiquitination K | N-glycosylation N | 8-task mean | 7-unit mean | Retained | Tasks < 0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ESM-2 | original | -0.0266 | +0.0149 | -0.1776 | -0.2488 | -0.2238 | -0.0663 | -0.0114 | -0.0666 | -0.1008 | -0.0814 | 1.000 | 7/8 |
| ESM-2 | relabel | -0.0247 | +0.0160 | -0.1715 | -0.2461 | -0.2245 | -0.0659 | -0.0060 | -0.0663 | -0.0986 | -0.0791 | 0.972 | 7/8 |
| ESM-2 | drop | -0.0265 | +0.0163 | -0.1869 | -0.2491 | -0.2241 | -0.0678 | -0.0105 | -0.0663 | -0.1019 | -0.0826 | 1.015 | 7/8 |
| ProtT5 | original | -0.0346 | +0.0133 | -0.2318 | -0.2414 | -0.2431 | -0.0735 | -0.0127 | -0.0631 | -0.1109 | -0.0921 | 1.000 | 7/8 |
| ProtT5 | relabel | -0.0326 | +0.0125 | -0.2185 | -0.2389 | -0.2448 | -0.0767 | -0.0058 | -0.0632 | -0.1085 | -0.0894 | 0.971 | 7/8 |
| ProtT5 | drop | -0.0348 | +0.0146 | -0.2429 | -0.2418 | -0.2435 | -0.0757 | -0.0116 | -0.0630 | -0.1124 | -0.0937 | 1.018 | 7/8 |
| interaction | original | -0.0152 | -0.0001 | -0.2616 | -0.2585 | -0.2817 | -0.1196 | -0.0294 | -0.0339 | -0.1250 | -0.1043 | 1.000 | 8/8 |
| interaction | relabel | -0.0155 | +0.0056 | -0.2452 | -0.2546 | -0.2828 | -0.1197 | -0.0261 | -0.0336 | -0.1215 | -0.1005 | 0.963 | 7/8 |
| interaction | drop | -0.0153 | +0.0021 | -0.2736 | -0.2586 | -0.2821 | -0.1229 | -0.0306 | -0.0336 | -0.1268 | -0.1063 | 1.020 | 7/8 |

### Step 2: ESM-2 by cell and source

| Cell | Source | Correction | Negatives flagged (3 splits) | 8-task mean | 7-unit mean | Retained | Tasks < 0 |
|---|---|---|---|---|---|---|---|
| TU | — | original | 0 | -0.1008 | -0.0814 | 1.000 | 7/8 |
| TU | union | relabel | 36828 | -0.0986 | -0.0791 | 0.972 | 7/8 |
| TU | union | drop | 36828 | -0.1019 | -0.0826 | 1.015 | 7/8 |
| TU | UniProt:primary | relabel | 951 | -0.1006 | -0.0811 | 0.997 | 7/8 |
| TU | UniProt:primary | drop | 951 | -0.1007 | -0.0814 | 0.999 | 7/8 |
| TU | UniProt:primary+similarity | relabel | 1582 | -0.1006 | -0.0811 | 0.996 | 7/8 |
| TU | UniProt:primary+similarity | drop | 1582 | -0.1008 | -0.0814 | 1.000 | 7/8 |
| TU | CPLM:primary | relabel | 17218 | -0.0990 | -0.0796 | 0.977 | 7/8 |
| TU | CPLM:primary | drop | 17218 | -0.1020 | -0.0828 | 1.018 | 7/8 |
| TU | PhosphoSitePlus:primary | relabel | 19366 | -0.1004 | -0.0810 | 0.995 | 7/8 |
| TU | PhosphoSitePlus:primary | drop | 19366 | -0.1006 | -0.0812 | 0.998 | 7/8 |
| TU | PhosphoSitePlus:low_throughput | relabel | 1009 | -0.1003 | -0.0809 | 0.994 | 7/8 |
| TU | PhosphoSitePlus:low_throughput | drop | 1009 | -0.1007 | -0.0813 | 0.999 | 7/8 |
| TT | — | original | 0 | +0.0437 | +0.0398 | 1.000 | 0/8 |
| TT | union | relabel | 24907 | +0.0389 | +0.0345 | 0.866 | 0/8 |
| TT | union | drop | 24907 | +0.0428 | +0.0387 | 0.973 | 0/8 |
| TT | UniProt:primary | relabel | 732 | +0.0436 | +0.0397 | 0.997 | 0/8 |
| TT | UniProt:primary | drop | 732 | +0.0436 | +0.0396 | 0.996 | 0/8 |
| TT | UniProt:primary+similarity | relabel | 1157 | +0.0434 | +0.0396 | 0.994 | 0/8 |
| TT | UniProt:primary+similarity | drop | 1157 | +0.0435 | +0.0396 | 0.994 | 0/8 |
| TT | CPLM:primary | relabel | 10252 | +0.0406 | +0.0364 | 0.914 | 0/8 |
| TT | CPLM:primary | drop | 10252 | +0.0427 | +0.0387 | 0.972 | 0/8 |
| TT | PhosphoSitePlus:primary | relabel | 14501 | +0.0422 | +0.0381 | 0.957 | 0/8 |
| TT | PhosphoSitePlus:primary | drop | 14501 | +0.0439 | +0.0400 | 1.005 | 0/8 |
| TT | PhosphoSitePlus:low_throughput | relabel | 779 | +0.0440 | +0.0401 | 1.007 | 0/8 |
| TT | PhosphoSitePlus:low_throughput | drop | 779 | +0.0438 | +0.0399 | 1.002 | 0/8 |
| UU | — | original | 0 | +0.0376 | +0.0384 | 1.000 | 0/8 |
| UU | union | relabel | 36828 | +0.0377 | +0.0385 | 1.002 | 0/8 |
| UU | union | drop | 36828 | +0.0388 | +0.0397 | 1.034 | 0/8 |
| UU | UniProt:primary | relabel | 951 | +0.0378 | +0.0386 | 1.006 | 0/8 |
| UU | UniProt:primary | drop | 951 | +0.0376 | +0.0385 | 1.002 | 0/8 |
| UU | UniProt:primary+similarity | relabel | 1582 | +0.0379 | +0.0388 | 1.011 | 0/8 |
| UU | UniProt:primary+similarity | drop | 1582 | +0.0377 | +0.0386 | 1.003 | 0/8 |
| UU | CPLM:primary | relabel | 17218 | +0.0388 | +0.0397 | 1.033 | 0/8 |
| UU | CPLM:primary | drop | 17218 | +0.0385 | +0.0395 | 1.027 | 0/8 |
| UU | PhosphoSitePlus:primary | relabel | 19366 | +0.0365 | +0.0371 | 0.967 | 0/8 |
| UU | PhosphoSitePlus:primary | drop | 19366 | +0.0378 | +0.0387 | 1.006 | 0/8 |
| UU | PhosphoSitePlus:low_throughput | relabel | 1009 | +0.0378 | +0.0387 | 1.007 | 0/8 |
| UU | PhosphoSitePlus:low_throughput | drop | 1009 | +0.0376 | +0.0385 | 1.002 | 0/8 |

### Enrichment: AUROC of the TU score for annotated vs other unrestricted test negatives (union; mean over three splits)

| Task | Baseline | ESM-2 | ESM-2 − baseline | ProtT5 | Interaction |
|---|---|---|---|---|---|
| Phosphorylation S/T | 0.636 | 0.616 | -0.0193 | 0.610 | 0.620 |
| Phosphorylation Y | 0.653 | 0.669 | +0.0156 | 0.663 | 0.663 |
| Acetylation K | 0.668 | 0.534 | -0.1331 | 0.509 | 0.493 |
| Methylation K/R | 0.506 | 0.409 | -0.0974 | 0.409 | 0.449 |
| Methylation R | 0.578 | 0.420 | -0.1574 | 0.178 | 0.244 |
| Sumoylation K | 0.614 | 0.568 | -0.0466 | 0.544 | 0.513 |
| Ubiquitination K | 0.590 | 0.597 | +0.0071 | 0.602 | 0.577 |
| N-glycosylation N | 0.925 | 0.860 | -0.0653 | 0.853 | 0.891 |
