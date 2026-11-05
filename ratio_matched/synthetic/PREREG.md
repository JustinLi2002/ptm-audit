# Synthetic ratio-matched control: fixed design (written 2026-11-01, before any run)

## Question
In the simulation, threshold training (T) changes both the donor entities and the achieved class ratio
(training positive rate ~0.56 at T = 10 against ~0.25 without a threshold). Does the reversal survive when
the training set without a threshold is given exactly T's positives and T's number of negatives?

## Fixed settings
- Generator, model and scoring: `analysis/synthetic_reversal.py` (generate, construct, fit_predict,
  measured_r2, homogeneity), unmodified, imported from the deposited repository.
- One point: threshold T = 10, r2_target = 0.60 (the highest fidelity of the published scan), 6,000 entities.
- Populations: seeds 0, 1, 2, 3, 4. Population, entity split (20% of entities held out) and the training
  construction of T and U follow the random-number sequence of `one_configuration()` exactly, so seeds 0-2
  of T and U must reproduce the published scan (`figures/figure6_sweep.tsv`, threshold 10 and 0, r2 0.60).
- Tests, identical for all three constructions within a population (rng seed + 99, as in the original):
  unrestricted test = construct(threshold 0) restricted to held-out entities (primary);
  threshold test = construct(threshold 10) restricted to held-out entities (secondary).

## Training constructions (training entities only; test sets untouched; no positive dropped)
- T: original rule, negatives only from entities with >= 10 positives, up to five per positive (pool-capped).
- U: original rule, negatives from all entities, five per positive.
- U_matched: all training positives (the same set as T and U) plus exactly as many negatives as T has,
  drawn uniformly without replacement from every unlabelled item of the training entities
  (rng seed + 1000). Positive count, negative count and positive rate therefore equal T's.

## Models
For each construction: base (covariate only), real (covariate + entity feature), perm (covariate +
permuted entity feature). fit_predict() as published: permutation mapping and initialisation seed are both
the population seed, hence shared across constructions within a population. Each model is fitted once and
scored on both tests. 3 constructions x 3 models x 5 populations = 45 fits.

## Outcomes
- Primary: d_real = AUROC(real) - AUROC(base of the same construction), unrestricted test.
- Secondary: the same for AUPRC; the same on the threshold test; real vs perm; change in
  rho(entity mean score, depth) on the unrestricted test.
- Recorded per run: positive and negative counts, positive rate, training entity count, share of training
  entities that are pure-positive and share of training items on them, absolute AUROC/AUPRC of all models
  on both tests, seeds, SHA-256 of both scripts. Predictions are saved.
- Populations are the replicates (n = 5); items and entities are not.

## Decision rule (applied in this order to the primary outcome)
1. Reproduction: T < 0 in all five populations and U > 0 in all five; otherwise inconclusive -> check
   implementation first.
2. U_matched > 0 in all five and mean(U_matched) closer to mean(U) than to mean(T) -> class ratio and
   sample size alone do not explain the contrast. If, in addition, U_matched - U has the same sign in all
   five populations, also state that matching changed the effect, and report effect sizes, not signs only.
3. mean(U_matched) <= 0, or mean(U_matched) closer to mean(T) than to mean(U) -> the mechanism attribution
   in the manuscript must be narrowed.
4. Anything else (U_matched of mixed sign) -> inconclusive.

## Limitation fixed in advance
U -> U_matched changes the class ratio, the training set size and the coverage of negatives at once, so no
"share due to ratio" can be computed. The directly controlled comparison is T against U_matched: equal
positives, equal counts, different donor entities. Conclusions are limited to this simulation point.
