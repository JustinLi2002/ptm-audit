# External-annotation check: results with UniProt and CPLM (PhosphoSitePlus pending)

Run on 2026-11-02 under PREREG.md, with no change to the pre-registered analysis. The status is **provisional**: two of
the three pre-registered sources are in. PREREG defines the primary external set as the union of all three, and
PhosphoSitePlus, which J.L. has to download, is not yet included. The missing source matters most for the two
phosphorylation tasks, where UniProt is currently the only source.

## Primary analysis (TU, ESM-2, relabel, union of sources): harm survives

| | Original | Corrected |
|---|---|---|
| 8-task mean ΔAUROC | −0.1008 | −0.0989 |
| 7-unit mean ΔAUROC | −0.0814 | −0.0794 |
| Retained share | — | 0.975 |
| Tasks below zero | 7/8 | 7/8 |

Under the decision rule (mean below zero, retained share ≥ 0.5, at least six of eight tasks below zero), the harm
survives. Step 3, retraining on cleaned negatives, is not triggered.

## Inputs

- Zenodo data record 21670043: rebuilt_datasets, replica_reconstruction, protein_disjoint_splits and
  crosseval_predictions_v2. All four MD5 checksums match the record.
- UniProtKB/Swiss-Prot release 2026_03, reference proteome UP000005640, reviewed entries. Fetched on 2026-11-02 through
  the REST search endpoint: 20,416 entries, with the modified-residue, glycosylation and cross-link features.
- CPLM Homo_sapiens.zip (16,739,171 bytes), downloaded on 2026-11-02. The file inside is dated 2021-06-13.
- Reference sequences: pdisjoint/proteins.fasta from the data record (16,336 proteins).

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
  - UniProt, experimental evidence: 42,366 sites matched. 497 were dropped for a residue mismatch, and 259 lie on
    proteins outside the reference.
  - CPLM: 220,617 sites matched. 443 were dropped for a residue mismatch, and 35,847 lie on proteins outside the
    reference.
- **Sequence drift between sources (diagnostic, post hoc).** For each CPLM site that matched on residue and has a
  full 15-residue window, the window in CPLM's own sequence was compared with ours. 215,816 windows are identical and
  60 differ, so sequence drift between the sources is negligible.

## What the results show

1. **The correction barely moves the harm.**
   - Across the three splits, 18,063 TU test negatives are externally annotated.
   - Relabelling them as positives moves the ESM-2 eight-task mean from −0.1008 to −0.0989.
   - Dropping them makes the mean slightly more negative (−0.1020).
   - ProtT5 and the interaction embedding behave the same way, with retained shares of 0.972 and 0.968.
   - Only one cell changes sign: the interaction embedding on phosphorylation Y, from −0.0001 to +0.0016 under
     relabelling.
2. **The external annotation is denser in the threshold-sampled negatives than in the unrestricted ones.**
   - In acetylation, sumoylation and ubiquitination, 8–16% of threshold-sampled negatives are externally annotated,
     against 4–10% of unrestricted negatives.
   - Within the unrestricted test negatives, those on donor proteins are annotated at a higher rate than those on
     non-donor proteins:
     - about 2.3 to 3.2 times in acetylation, sumoylation, ubiquitination and N-glycosylation;
     - more in the two methylation tasks, at low absolute rates;
     - 1.8 times in phosphorylation S/T.
   - So annotated non-sites are recorded at higher rates on well-studied proteins, which are the proteins the threshold selects. The Methods
     sentence that the residual noise "need not be equal" between constructions can now be given a direction and a
     size.
   - Phosphorylation Y is the exception: there the rate is higher off donor proteins (0.70% against 0.38%).
3. **The direct form of the alternative explanation runs the wrong way.**
   - If the feature pushed true but unannotated sites up, the annotated negatives would rank higher under ESM-2 than
     under the baseline.
   - In six of eight tasks they rank lower. For example, the enrichment AUROC falls from 0.668 to 0.534 in
     acetylation K and from 0.578 to 0.420 in methylation R.
   - ESM-2 ranks them higher only in phosphorylation Y (+0.049), the one task where TU ESM-2 is not harmful, and
     marginally in ubiquitination K (+0.007).
   - This pattern is consistent with the shortcut described in the paper: a threshold-trained feature learns that
     negatives lie on donor proteins, scores negatives on well-studied proteins low, and the annotated negatives
     are recorded at higher rates on exactly those proteins. This reading is an interpretation, not a test.
4. **The in-construction cells are hardly affected.**
   - TT: relabelling lowers the ESM-2 gain from +0.0437 to +0.0405 (retained share 0.911). A small part of the gain on
     threshold-sampled test sets comes from ranking annotated negatives low.
   - UU: the gain is unchanged or slightly larger.

## Deviations from PREREG.md and additions to it

- **Python flags.** Python was run with `-E -P` rather than `-I`, because numpy is installed in the user site-packages.
  With these flags, neither the script directory nor the working directory is on the import path. This has no effect
  on the results.
- **Union.** The union currently contains UniProt and CPLM only.
- **Additions (descriptive, not pre-registered):**
  - validation-partition counts and the annotated share of the positives in step1_counts.tsv;
  - the CPLM window diagnostic;
  - UniProt similarity-only evidence as a separate source set;
  - ProtT5 and the interaction embedding in the enrichment analysis, which PREREG specified for the baseline and
    ESM-2 only.
- **CPLM evidence.** CPLM evidence includes "Dat." (integrated from other databases). As PREREG specifies, all sites
  were used. No experimental-only CPLM tier was run; one would be labelled post hoc.

## Pending

PhosphoSitePlus. J.L. downloads these site datasets into `data/psp/`:
- Phosphorylation_site_dataset;
- Acetylation_site_dataset;
- Methylation_site_dataset;
- Ubiquitination_site_dataset;
- Sumoylation_site_dataset.

Then rerun build_sites.py and extannot_check.py (commands in scripts/README.txt). The human rows are read
automatically, the primary set becomes the three-source union, and the low-throughput tier (LT_LIT ≥ 1) is added
automatically.

## Files in this folder

- `report.txt`: checks, matching counts and the primary verdict.
- `sites_report.txt`, `sites.tsv`: the parsed external sites.
- `step1_counts.tsv`: counts by task, split, construction, partition and source.
- `step1_test_negatives_by_donor_status.tsv`: unrestricted test negatives on donor and non-donor proteins.
- `step2_rescoring_by_split.tsv`, `step2_summary.tsv`: the rescoring, per split and summarised.
- `enrichment_by_split.tsv`, `enrichment_summary.tsv`: the enrichment analysis, per split and summarised.

## Tables (generated from the TSV files by make_results_md.py)

### Step 1: share of negatives that are externally annotated (union; mean over three splits)

| Task | T train neg % | T test neg % | U train neg % | U test neg % | test pos % | U test neg % (donor proteins) | U test neg % (non-donor proteins) |
|---|---|---|---|---|---|---|---|
| Phosphorylation S/T | 0.31 | 0.32 | 0.25 | 0.25 | 14.60 | 0.31 | 0.17 |
| Phosphorylation Y | 0.25 | 0.28 | 0.39 | 0.45 | 7.10 | 0.38 | 0.70 |
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
| ESM-2 | relabel | -0.0264 | +0.0158 | -0.1715 | -0.2461 | -0.2245 | -0.0659 | -0.0060 | -0.0663 | -0.0989 | -0.0794 | 0.975 | 7/8 |
| ESM-2 | drop | -0.0266 | +0.0152 | -0.1869 | -0.2491 | -0.2241 | -0.0678 | -0.0105 | -0.0663 | -0.1020 | -0.0828 | 1.017 | 7/8 |
| ProtT5 | original | -0.0346 | +0.0133 | -0.2318 | -0.2414 | -0.2431 | -0.0735 | -0.0127 | -0.0631 | -0.1109 | -0.0921 | 1.000 | 7/8 |
| ProtT5 | relabel | -0.0342 | +0.0137 | -0.2185 | -0.2389 | -0.2448 | -0.0767 | -0.0058 | -0.0632 | -0.1085 | -0.0895 | 0.972 | 7/8 |
| ProtT5 | drop | -0.0345 | +0.0134 | -0.2429 | -0.2418 | -0.2435 | -0.0757 | -0.0116 | -0.0630 | -0.1125 | -0.0939 | 1.019 | 7/8 |
| interaction | original | -0.0152 | -0.0001 | -0.2616 | -0.2585 | -0.2817 | -0.1196 | -0.0294 | -0.0339 | -0.1250 | -0.1043 | 1.000 | 8/8 |
| interaction | relabel | -0.0150 | +0.0016 | -0.2452 | -0.2546 | -0.2828 | -0.1197 | -0.0261 | -0.0336 | -0.1219 | -0.1010 | 0.968 | 7/8 |
| interaction | drop | -0.0152 | +0.0004 | -0.2736 | -0.2586 | -0.2821 | -0.1229 | -0.0306 | -0.0336 | -0.1270 | -0.1066 | 1.022 | 7/8 |

### Step 2: ESM-2 by cell and source

| Cell | Source | Correction | Negatives flagged (3 splits) | 8-task mean | 7-unit mean | Retained | Tasks < 0 |
|---|---|---|---|---|---|---|---|
| TU | — | original | 0 | -0.1008 | -0.0814 | 1.000 | 7/8 |
| TU | union | relabel | 18063 | -0.0989 | -0.0794 | 0.975 | 7/8 |
| TU | union | drop | 18063 | -0.1020 | -0.0828 | 1.017 | 7/8 |
| TU | UniProt:primary | relabel | 951 | -0.1006 | -0.0811 | 0.997 | 7/8 |
| TU | UniProt:primary | drop | 951 | -0.1007 | -0.0814 | 0.999 | 7/8 |
| TU | UniProt:primary+similarity | relabel | 1582 | -0.1006 | -0.0811 | 0.996 | 7/8 |
| TU | UniProt:primary+similarity | drop | 1582 | -0.1008 | -0.0814 | 1.000 | 7/8 |
| TU | CPLM:primary | relabel | 17218 | -0.0990 | -0.0796 | 0.977 | 7/8 |
| TU | CPLM:primary | drop | 17218 | -0.1020 | -0.0828 | 1.018 | 7/8 |
| TT | — | original | 0 | +0.0437 | +0.0398 | 1.000 | 0/8 |
| TT | union | relabel | 10903 | +0.0405 | +0.0363 | 0.911 | 0/8 |
| TT | union | drop | 10903 | +0.0426 | +0.0385 | 0.968 | 0/8 |
| TT | UniProt:primary | relabel | 732 | +0.0436 | +0.0397 | 0.997 | 0/8 |
| TT | UniProt:primary | drop | 732 | +0.0436 | +0.0396 | 0.996 | 0/8 |
| TT | UniProt:primary+similarity | relabel | 1157 | +0.0434 | +0.0396 | 0.994 | 0/8 |
| TT | UniProt:primary+similarity | drop | 1157 | +0.0435 | +0.0396 | 0.994 | 0/8 |
| TT | CPLM:primary | relabel | 10252 | +0.0406 | +0.0364 | 0.914 | 0/8 |
| TT | CPLM:primary | drop | 10252 | +0.0427 | +0.0387 | 0.972 | 0/8 |
| UU | — | original | 0 | +0.0376 | +0.0384 | 1.000 | 0/8 |
| UU | union | relabel | 18063 | +0.0389 | +0.0399 | 1.038 | 0/8 |
| UU | union | drop | 18063 | +0.0386 | +0.0395 | 1.029 | 0/8 |
| UU | UniProt:primary | relabel | 951 | +0.0378 | +0.0386 | 1.006 | 0/8 |
| UU | UniProt:primary | drop | 951 | +0.0376 | +0.0385 | 1.002 | 0/8 |
| UU | UniProt:primary+similarity | relabel | 1582 | +0.0379 | +0.0388 | 1.011 | 0/8 |
| UU | UniProt:primary+similarity | drop | 1582 | +0.0377 | +0.0386 | 1.003 | 0/8 |
| UU | CPLM:primary | relabel | 17218 | +0.0388 | +0.0397 | 1.033 | 0/8 |
| UU | CPLM:primary | drop | 17218 | +0.0385 | +0.0395 | 1.027 | 0/8 |

### Enrichment: AUROC of the TU score for annotated vs other unrestricted test negatives (union; mean over three splits)

| Task | Baseline | ESM-2 | ESM-2 − baseline | ProtT5 | Interaction |
|---|---|---|---|---|---|
| Phosphorylation S/T | 0.693 | 0.680 | -0.0134 | 0.688 | 0.689 |
| Phosphorylation Y | 0.666 | 0.715 | +0.0487 | 0.701 | 0.719 |
| Acetylation K | 0.668 | 0.534 | -0.1331 | 0.509 | 0.493 |
| Methylation K/R | 0.506 | 0.409 | -0.0974 | 0.409 | 0.449 |
| Methylation R | 0.578 | 0.420 | -0.1574 | 0.178 | 0.244 |
| Sumoylation K | 0.614 | 0.568 | -0.0466 | 0.544 | 0.513 |
| Ubiquitination K | 0.590 | 0.597 | +0.0071 | 0.602 | 0.577 |
| N-glycosylation N | 0.925 | 0.860 | -0.0653 | 0.853 | 0.891 |
