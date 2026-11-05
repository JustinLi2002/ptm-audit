# Real-data ratio-matched control: fixed design (written 2026-11-01, before any run)

Follows the synthetic control (../ratio_matched/), whose outcome fell in row 2 of its decision rule
("ratio and sample size alone do not explain the contrast"), for which the plan's next step is validation on
real tasks.

## Tasks and feature
- Acetylation K and phosphorylation S/T. Both are pool-limited under the threshold (published training
  positive rates 0.43 and 0.27, against 0.167 in the unrestricted reconstruction), so the matched
  construction differs from the unrestricted one in class ratio and size. Methylation K/R is not run.
- ESM-2 only (frozen, mean-pooled, as in the manuscript).

## Constructions
- T (threshold-sampled, `replica/`) and U (unrestricted, `rebuilt/`): existing runs in `pdisjoint_runs_v2/`,
  reused unchanged (same script, partitions, initialisations and permutation seed).
- Matched (`matched_v4/split{s}/`): for each task, partition seed s in {0, 1, 2} and each of the train and val
  partitions, all positives of the unrestricted reconstruction in that partition plus a uniform random
  subsample, without replacement, of its unrestricted negatives, of exactly the threshold-sampled negative count
  in that partition (pandas sample, random_state 1000 + s for train and 2000 + s for val). If the two
  constructions' positive sets in a partition differ, the count is scaled to keep T's negative-to-positive
  ratio, and this is reported. A uniform subsample of the uniformly sampled unrestricted negatives is a uniform
  sample of the unrestricted candidate pool. Test partitions are untouched: every model is scored on both
  existing test constructions.
- Validation (early stopping) uses the matched val partition, as T and U each use their own.

## Models
`ptm-audit/train_pdisjoint.py`, run through `train_matched.py`, which changes only the training-data path and
the output directories. Conditions baseline, ESM-2 real ("ppi" with --feat esm) and ESM-2 permuted
("shuffled", permutation seed 42 as in the existing runs); two initialisations (0, 1) per partition; 3
conditions x 3 partitions x 2 initialisations = 18 fits per task, 36 in all. The loss is weighted by the
training class ratio, as in every run of the paper.

## Outcomes
- Primary: ESM-2 dAUROC = AUROC(real) - AUROC(baseline of the same construction) on the unrestricted test,
  ensemble of the two initialisations, per partition and averaged over the three partitions, per task.
- Secondary: dAUPRC; the same on the threshold-sampled test; real vs permuted; training positive rate,
  negative count and pure-positive share of each construction.
- Replicates are the three partitions within a task; tasks are reported separately and not pooled into a test.

## Decision rule, per task, in this order
1. If the existing T is not negative on average, or the existing U not positive, the task cannot address the
   question; report it descriptively.
2. Matched > 0 in all three partitions and mean(matched) closer to mean(U) than to mean(T): the class ratio
   and training-set size alone do not explain the contrast in this task.
3. mean(matched) <= 0, or mean(matched) closer to mean(T) than to mean(U): the mechanism attribution must be
   narrowed for this task.
4. Anything else: mixed; report effect sizes.

## Limitation fixed in advance
U -> matched changes the class ratio, the training-set size and the coverage of negatives together, and the
loss weighting already compensates part of the class ratio, so no share of the effect is attributed to the
ratio. The directly controlled comparison is T against matched: equal positives, equal negative counts,
different donor proteins. Conclusions are limited to these two tasks and this feature.
