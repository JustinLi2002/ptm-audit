# Independent-source replication, mouse phosphorylation S/T: fixed design

FIXED 2026-11-02: approved by Justin as written, before any of the following was done:
- reference sequences mapped;
- any construction statistic beyond the preflight counts computed;
- any model trained.

Changes after this point are reported as deviations.

## Question

Does the same negative-sampling intervention change the measured contribution of ESM-2 in the same way when applied to phosphorylation data from an independent experimental source? Here "the same intervention" means drawing negatives only from proteins with at least ten sites, against drawing them from all site-bearing proteins.

This is an intervention we apply to the data. It is not a claim that the original study used a threshold or made any error.

## Data

**Source.** Giansanti et al., Nat Methods 2022;19:803–11, doi:10.1038/s41592-022-01526-y; PRIDE PXD030983. We use Supplementary Table 2, sheet "Supplementary Table 2 - Organs": 41 tissues of C57BL/6N mice, MaxQuant 1.6.0.13 phosphosite table.
- Workbook sha256: 9f397a5a68310fcbcc0628729d97ac2b701d5de37e8dee7a41c33d98671e0e08.
- The mPDAC cell-line sheet is not used.
- The source is independent of DeepMVP's PTMAtlas in species, laboratory, samples and processing. It shares the feature model (ESM-2) and our training code.

**Site filters**, as in the preflight:
- localization probability ≥ 0.75;
- an intensity in at least one of the 41 tissues;
- a single leading protein;
- the centre of the 31-residue sequence window equals the reported residue.

Duplicates are removed on accession, position and residue. The preflight gives:
- 29,039–29,054 S/T sites (the two counts differ by the deduplication detail; the mapped count below is final);
- 564 proteins with merged S/T/Y depth ≥ 10.

**Reference sequences.** UniProtKB/Swiss-Prot *Mus musculus*, canonical sequences only, as for the human reconstruction. The release and download date are recorded.
- The original search used UP000000589 (March 2018), so positions are not reused. Each site is placed by matching its 31-residue window, with "_" padding at the termini treated as sequence ends, exactly and uniquely to the Swiss-Prot canonical entry of an accession in the row's "Proteins" field. The leading accession's own entry is preferred.
- Sites without a unique exact match are dropped.
- If fewer than 80% of S/T sites map, we stop and report before building anything.

## Constructions

`rebuild_datasets.py` logic, unchanged except for inputs.
- **Positives:** all mapped S/T sites.
- **Depth:** mapped S/T/Y sites per protein. Phosphorylation S/T and Y are merged, as in the human threshold group.
- **T (threshold-sampled):** candidate negatives are all S/T residues that are not positives, on proteins with depth ≥ 10.
- **U (unrestricted):** candidate negatives are the same, on all proteins with at least one S/T positive.
- **Sampling:** uniform, 5 negatives per positive, capped by the pool, seed 0. These are the human settings; none is tuned here.
- **Matched:** as in `restructure_20261101/v4/ratio_matched_real/PREREG_real.md`. In each training and validation partition it has U's positives plus a uniform subsample of U's negatives, of T's negative count in that partition (random_state 1000 + s for training, 2000 + s for validation).
- **Splits:**
  - CD-HIT -c 0.4 -n 2 over all proteins in U;
  - `make_pdisjoint_split.py` with test 0.10, validation 0.10, seeds 0, 1, 2;
  - T, U and matched share the protein-to-partition assignment.

Recorded before training, as the construction-level prediction: the pure-positive share of training sites in T, U and matched. Pure-positive means sites on proteins without a negative. Human reference values: phosphorylation S/T T 3.9%; acetylation K T 20.6%, matched 2.2%, U 0.05%. Across human tasks the harm in TU grew with this share, so a T share well above 3.9% predicts a larger TU harm than in human phosphorylation S/T.

## Models

**Feature.** ESM-2 650M (esm2_t33_650M_UR50D), frozen, mean-pooled over the full sequence (`extract_esm.py`), uncentred.

**Training.** `train_pdisjoint.py` with its current hyperparameters, run through a wrapper that changes only the data paths.
- Conditions: sequence-only baseline; ESM-2 real; ESM-2 permuted across proteins (seed 42).
- Two initialisations per split.
- 3 training constructions × 3 splits × 3 conditions × 2 initialisations = 54 fits.
- Every model is scored on the T test and on the U test of its split. Neither test is used for any choice.
- On HRP, GPU jobs exclude the Tesla M10 nodes, and failures fail the job.

## Outcomes

**Primary statistic:** ESM-2 ΔAUROC (real minus baseline), using the ensemble of the two initialisations, averaged over three splits.
- 95% interval from the protein bootstrap of `bootstrap_20261102/crosseval_bootstrap_v2.py`: joint resampling across splits, the split mean as the statistic, 2,000 resamples.
- An interval is *above* or *below* zero when it excludes zero on that side.

**Decision rule.** The first row that applies is reported, with its label, whichever it is.

| Row | Condition | Label |
|---|---|---|
| A | TT above zero and TU below zero | Direct replication: the feature helps on the threshold-sampled test and harms on the unrestricted test |
| B | TU below zero and UU above zero, and A does not hold | Training-construction effect replicated; the threshold-sampled test does not hide it |
| C | TU below zero, and neither A nor B holds | Harm on the unrestricted test after threshold training, without the reversal |
| D | TU not below zero | No harm: not replicated |

Reported alongside, without changing the row:
- per-split signs in each cell;
- ΔAUPRC;
- real against permuted ESM-2 (split pairs, both tests);
- the construction-level prediction above.

**Secondary (matched).** ESM-2 ΔAUROC on the U test after matched training, against T and U training.

Read as in the real-data control:
- matched close to U and different from T: class ratio and size alone do not account for the T result;
- matched close to T: they may;
- otherwise: unresolved.

## Fixed in advance

No change after any model result:
- the threshold (10);
- the ratio (5);
- the filters;
- the window;
- the reference;
- the split fractions;
- the feature;
- the training settings.

No other threshold or task is tried. A technically failed run is rerun with the same configuration and reported.

## Not claimed

The negatives are unobserved sites, not sites shown to be unmodified.

A second step, restricting candidate negatives to residues covered by observed peptides, needs the search-level outputs. It would get its own fixed design. Its negatives would be called candidate negatives supported by observation.
