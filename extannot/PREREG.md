# External-annotation check of the unrestricted test negatives (pre-registration)

Written 2026-11-02, before any external database file was downloaded or any hit counted. Every analysis below is fixed
here; anything added after the first result is seen will be labelled post hoc.

## Question

The harm of protein-level features under threshold-sampled training (TU: ESM-2 mean ΔAUROC −0.101 over eight tasks) is
measured against unrestricted test negatives, which are unannotated candidate residues, not observed non-sites. The
alternative explanation is that some of these negatives are true sites recorded in other databases. Under this
explanation, the threshold-trained model with the feature would be scoring them highly and being penalised for being
right. The original DeepMVP protocol removed sites annotated in PhosphoSitePlus, UniProt and PLMD from the negative pool.
Our reconstructions did not do this (Methods), and because the two constructions draw negatives from different proteins,
the noise it leaves need not be equal.

This check uses external annotation rather than the depth-based missingness model of Supplementary Note 6.

## Sources and evidence filters (human only)

| Source | Version | Primary filter | Sensitivity tier |
|---|---|---|---|
| PhosphoSitePlus | Files downloaded by J.L. (date recorded at download) | Any site in the organism = human rows | Low-throughput only (LT_LIT ≥ 1) |
| UniProtKB/Swiss-Prot | Release 2026_03 (2 Sep 2026), reference proteome UP000005640, reviewed entries | Features with experimental evidence: ECO:0000269 or ECO:0007744 | Adds by-similarity (ECO:0000250). Sequence-analysis-only evidence (ECO:0000255) is never used, because it would flag every N-glycosylation sequon |
| CPLM (the CUCKOO group's successor to PLMD) | Homo sapiens archives, site last updated 10 Feb 2025 | All sites of the mapped types | — |

The union of the three sources, each under its primary filter, is the primary external set.

## Task → modification types

| Task | PhosphoSitePlus | UniProt feature | CPLM |
|---|---|---|---|
| Phosphorylation S/T | Phosphorylation, S/T | Modified residue: phosphoserine, phosphothreonine | — |
| Phosphorylation Y | Phosphorylation, Y | Modified residue: phosphotyrosine | — |
| Acetylation K | Acetylation, K | Modified residue: N6-acetyllysine | Acetylation |
| Methylation K/R (released as methylation K; contains K and R sites) | Methylation, K and R | Modified residue: methyl-, dimethyl-, trimethyllysine and methylated arginines | Methylation (K) |
| Methylation R | Methylation, R | Modified residue: methylated arginines (omega-N-methyl, symmetric and asymmetric dimethyl) | — |
| Ubiquitination K | Ubiquitination | Cross-link: Glycyl lysine isopeptide, interchain with G-Cter in ubiquitin | Ubiquitination |
| Sumoylation K | Sumoylation | Cross-link: Glycyl lysine isopeptide, interchain with G-Cter in SUMO1/2/3 | Sumoylation |
| N-glycosylation N | — | Glycosylation: N-linked (GlcNAc...) asparagine | — |

Task file names are checked against this table before any counting.

## Matching

- A database site matches a candidate site when all of these hold:
  - same UniProt accession;
  - same position;
  - the residue at that position in our reference sequence is the task residue, and it agrees with the residue the database states.
- Isoform-specific accessions (with a suffix other than the canonical one) are not used.
- Sites whose residue does not match are dropped and counted per source.
- Sites already labelled positive in our sets are unaffected; only negatives can be hit.

## Step 1: counts (no models)

For each task, split, construction (T, U) and partition (train, test), we count the externally annotated negatives and
their share of all negatives. Counts are given by source and for the union.

Unrestricted test negatives are also split by donor status:
- on proteins that are donors under the threshold (at least ten merged-type sites);
- on proteins that are not, i.e. the negatives that only the unrestricted construction contains.

## Step 2: rescoring with the deposited predictions (no retraining)

The input is the deposited per-site predictions of the threshold-trained ensembles on the unrestricted test (TU), for:
- the sequence-only baseline;
- ESM-2 (the primary feature);
- the interaction embedding and ProtT5 (secondary features).

Two corrections are applied to the externally annotated negatives:
- **relabel** (primary): they become positives;
- **drop** (secondary): they are removed from the test set.

The quantity reported is ΔAUROC (feature minus baseline). It is given per task as the mean over three splits, as the mean
over eight tasks, and as the mean over seven units, which collapses the two methylation tasks as in Supplementary Note 6.
For each result we also report how many tasks are negative, and the retained share (corrected seven-unit mean divided by
the original).

TT and UU are rescored the same way for reference.

### Enrichment among the unrestricted test negatives (secondary)

Within each task and split, AUROC of the TU score for separating externally annotated negatives from the other
negatives, computed for the baseline and for ESM-2. If this AUROC is higher for ESM-2 than for the baseline, the
negatives that the feature pushes up are disproportionately known sites, which is the alternative explanation in its
direct form.

## Decision rule (primary analysis: relabel, union of sources, ESM-2, TU)

- **Harm survives:** the eight-task mean stays below zero, the retained share is at least 0.5, and at least six of the
  eight task means stay below zero.
- **Partly explained:** the mean stays below zero, but the retained share is below 0.5 or fewer than six tasks stay
  negative.
- **Explained:** the mean is zero or above.

Step 3, retraining on cleaned training negatives, is run only if the result is "partly explained" or "explained", and only
after discussion.

## Limits stated in advance

- The external databases are themselves deeper for well-studied proteins. This check can find missing positives only
  where someone has looked, so it is weaker exactly where Supplementary Note 6's concern is largest.
- PhosphoSitePlus includes high-throughput mass-spectrometry sites of mixed reliability, which is why the low-throughput
  tier is reported.
- Relabelling treats every database site as a true site in the context of the benchmark; this is the adverse
  assumption.
- No interval estimates are pre-specified; signs, means and counts are reported as above.
