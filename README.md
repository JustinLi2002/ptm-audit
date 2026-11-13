# ptm-audit

Diagnostics for label leakage through protein-level feature channels in PTM site
prediction. Companion code for "Negative sampling determines the sign of protein-level feature contributions in
modification site prediction" (Li & Yu, manuscript in preparation).

## Current manuscript and what v1.6 adds

The manuscript was restructured after v1.5. The sections below this one were written for v1.5: they use that draft's
numbering and, in places, stronger wording than the manuscript now uses ("channel", "key into a lookup table", "natural
negatives"). Where they differ, the manuscript and this section take precedence.

Names used in files: T, threshold-sampled (`replica`); U, unrestricted (`rebuilt`); M, ratio-matched (`matched`). A cell
is the training construction followed by the test construction, so TU is trained on T and tested on U. A split is one of
the three protein-disjoint replicates; a partition is its training, validation or test part.

### Display items

| Item | Script | Reads or writes |
|---|---|---|
| Fig. 1 | `figures/v2026_11/fig1_construction_v8.py` | `results/unobserved_weight.txt` |
| Fig. 2 | `figures/v2026_11/fig2_reversal_v7.py` | `results/crosseval_verify.tsv`; panels g,h `ratio_matched/real/runs_real.tsv` |
| Fig. 3 | `figures/v2026_11/fig3_explanation_v4.py` | `results/section3_20260920/`, `results/crosseval_verify.tsv` |
| Fig. 4 | `figures/v2026_11/fig4_synthetic_v5.py` | `figures/figure6_sweep.tsv`; panel c `ratio_matched/synthetic/runs.tsv` |
| Table 1 | `analysis/table1_cells.py` | a: `results/crosseval_verify.tsv`, `results/crosseval_ci_v2.tsv`; b: `mouse/results/mouse_results.tsv` |
| Supp. Fig. 1 | `analysis/plot_corrected.py` | `results/section3_20260920/figureS4_corrected.png` |
| Supp. Fig. 2 | `analysis/plot_perm_margin.py` | `figures/figure_S5.png` |
| Supp. Fig. 3 | `figures/v2026_11/suppfig3_dose_response.py` | `results/crosseval_verify.tsv` |
| Supp. Fig. 4 | `figures/v2026_11/suppfig4_variance.py` | `results/icc_audit.tsv`, `figures/v2026_11/data_hpc_20261101/` |
| Supp. Fig. 5 | `analysis/make_figure_s2.py --results results/section3_20260920` | `results/section3_20260920/figure_s2_fn_sensitivity.png` |
| Supp. Fig. 6 | `make_figures.py` | `figures/figure2.png` |
| Supp. Fig. 7 | `make_figures.py` | `figures/figureS1.png` |
| Supp. Table 1 | `audit_ptm_benchmark.py` section T3, `analysis/benchmark_homogeneity.py` | |
| Supp. Tables 2-4 | `analysis/crosseval_summary.py [--all]`, `analysis/crosseval_verify.py` | `results/crosseval_verify.tsv` |
| Supp. Table 5 | `analysis/section3_corrected.py` | `results/section3_20260920/ridge_summary.tsv` |
| Supp. Table 6 | `train_alldata.py`, `analysis/summarize_alldata.py`, `analysis/s9_auprc_v2.py` | |
| Supp. Table 7 | `analysis/fn_sensitivity.py`, `analysis/recompute_methyl_sensitivity.py`, `analysis/assemble_section3_sensitivity.py`, `analysis/make_table_s10.py` | `results/section3_20260920/` |
| Supp. Table 8 | `analysis/benchmark_homogeneity.py` | `results/cross_benchmark.tsv` |
| Supp. Tables 9, 12, 13 | `audit_ptm_benchmark.py`, `restricted_eval.py` | |
| Supp. Tables 10, 11 | `analysis/cluster_bootstrap.py` | `results/cluster_bootstrap_out.txt` |
| Supp. Table 14 | `ratio_matched/real/analyse_real.py` | `ratio_matched/real/runs_real.tsv`, `ratio_matched/real/matched_v4/matched_summary.tsv` |
| Supp. Table 15 | `extannot/scripts/extannot_check.py` | `extannot/results/` |
| Supp. Table 16 | `analysis/crosseval_bootstrap_v2.py` | `results/crosseval_ci_v2.tsv` |
| Supp. Table 17 | `mouse/scripts/analyse_mouse.py` | `mouse/results/` |

Supplementary Note 4 (ratio-matched simulation) uses `ratio_matched/synthetic/`, Note 6 (external annotation) `extannot/`,
and Note 9 (mouse) `mouse/`.

### Added in v1.6

| Path | What it is |
|---|---|
| `analysis/crosseval_bootstrap_v2.py`, `results/crosseval_ci_v2.*`, `docs/crosseval_bootstrap_v2.txt` | The interval for the reported statistic, the mean over splits of the per-split difference. It supersedes `analysis/crosseval_bootstrap.py` and `results/crosseval_ci.tsv`, which pool the splits and are kept for provenance |
| `analysis/table1_cells.py`, `results/table1_cells.tsv` | Table 1 from the deposited tables |
| `ratio_matched/synthetic/` | Ratio-matched control in the simulation, with its PREREG |
| `ratio_matched/real/` | Ratio-matched control on acetylation K and phosphorylation S/T, with its PREREG |
| `extannot/` | External-annotation check, aggregate tables only; see `extannot/README.txt` for access to the databases |
| `mouse/` | Independent mouse dataset; see `mouse/README.md` |
| `figures/v2026_11/` | Scripts and outputs of Figs 1-4 and Supplementary Figs 3-4 |

The human ratio-matched training sets are not archived as files because they are reproducible: with the three folders
`replica/`, `rebuilt/` and `pdisjoint/` from the data record under `PTM_AUDIT_BASE`, `ratio_matched/real/build_matched.py`
gives the six files in `ratio_matched/real/matched_v4/manifest.sha256`, byte for byte.

### Changes made when copying into this repository

The scripts ran from other folders. Compared with the versions that ran, only the following differ:
- `.py`, `.sh` and `.slurm` files have LF line endings;
- hard-coded local paths became paths relative to the repository root, in `ratio_matched/synthetic/{ratio_matched,summarise}.py`, `ratio_matched/real/{analyse_real,train_matched}.py`, `mouse/scripts/analyse_mouse.py`, `mouse/scripts/hrp/train_mouse.py` and `figures/v2026_11/{figstyle,fig2_reversal_v7,fig4_synthetic_v4,fig4_synthetic_v5}.py`;
- `ratio_matched/real/analyse_real.py` takes its input folder as an optional argument;
- `ratio_matched/real/build_matched.py` reads `PTM_AUDIT_BASE` and always writes LF line endings;
- the login host in the fetch scripts is a placeholder, `USER@HPC_LOGIN_HOST`.

The SHA-256 of `ratio_matched.py` recorded in `ratio_matched/synthetic/runs.tsv` is that of the version that ran. The
slurm files keep the cluster paths they ran with, as in `slurm/`. `PREREG_mouse.md` and all tables are byte-identical to
the originals.

Checked with Python 3.14.3, numpy 2.5.1, pandas 3.0.3, scipy 1.18.0, scikit-learn 1.9.0 and matplotlib 3.11.0: from these
paths, the six figures above are pixel-identical to the manuscript's, `analyse_real.py` and `summarise.py` reproduce
their recorded outputs byte for byte, and `analyse_mouse.py`, run on the bundles in the data record, reproduces
`mouse/results/mouse_results.tsv` and `summary.txt` byte for byte. The figure scripts import from their own folder, so run
them as `python figures/v2026_11/<script>.py`.

## What the leakage is

Negative samples in some PTM benchmarks are drawn only from proteins carrying at
least a threshold number of annotated sites. Proteins below the threshold
contribute positives and no negatives, so every one of their sites carries the
same label. Any feature that is constant within a protein — an interaction
embedding, a kinase-substrate prior, a frozen protein language model embedding —
can then act as a key into a lookup table, provided the same protein appears in
both the training and test partitions.

Two conditions must therefore hold for the channel to be exploitable: a
permissive sampling rule, and a split that places the same protein on both
sides. Correcting either one closes it. Correcting only the split still leaves
the training set distorted.

The distortion has a direction. The threshold makes annotation depth almost
determine the label — proteins below it are all-positive — and a protein-level
embedding partially encodes annotation depth, because both interaction degree
and modification counts scale with how heavily a protein has been studied. A
model trained this way learns "shallowly annotated → positive", which is the
reverse of the real relationship. On unseen proteins under a natural negative
distribution the channel therefore ranks systematically backwards, and the real
embedding does *more* damage than a permuted vector of the same shape.

The convention is a choice rather than a necessity, and the field is split on
it. MusiteDeep draws negatives from every candidate residue of every sequence it
holds, and 0.0% of its 3.0 million training sites lie on label-homogeneous
proteins, against 4.4–40.8% in the release analysed here
(`analysis/benchmark_homogeneity.py`, Supplementary Table S12).

The remedy is not to abandon the threshold. Drawing negatives from well-assayed
proteins is defensible for assembling a *training* set, since a residue in a
protein that was never assayed is not evidence of absence. What does not follow
is that the same rule should assemble the partition on which the model is
judged.

## Figure and table map (draft numbering, up to v1.5)

Output filenames carry the numbering the figures had while the manuscript was
being drafted, which runs one ahead of the final numbering for Figures 3 to 5.
Renaming them would break the `savefig` paths in three scripts, so the mapping is
recorded here instead.

| Manuscript | File | Produced by |
|---|---|---|
| Figure 1 | `figures/figure1.png` | hand-drawn |
| Figure 2, annotation depth | `results/section3_20260920/figure2_corrected.png` | `analysis/plot_corrected.py`, from tables written by `analysis/section3_corrected.py` |
| Figure 3, cross-evaluation | `figures/figure4.png` | `analysis/make_figure_merged.py` |
| Figure 4, variance structure | `figures/figure5.png` | `analysis/make_figures_v2.py`, `variance_figure` |
| Figure 5, synthetic sweep | `figures/figure6.png` | `analysis/synthetic_reversal.py --mode sweep` |
| Supplementary Figure S1, permutation control | `figures/figureS1.png` | `make_figures.py` |
| Supplementary Figure S2, false-negative sensitivity | `results/section3_20260920/figure_s2_fn_sensitivity.png` | `analysis/make_figure_s2.py --results results/section3_20260920` |
| Supplementary Figure S3, identity baseline | `figures/figure2.png` | `make_figures.py` |
| Supplementary Figure S4, depth by candidate count | `results/section3_20260920/figureS4_corrected.png` | `analysis/plot_corrected.py` |
| Supplementary Figure S5, permutation margin | `figures/figure_S5.png` | `analysis/plot_perm_margin.py` |

Supplementary Figure S3 was Figure 2 of an earlier draft and moved to the
supplement when the Results were reordered; its file name did not follow.
Figures 1, S1, S2 and S3 also have `.pdf` companions in the same directory.

The Section 3 revision (2026-09-20) replaced Figure 2 and Supplementary Figures
S2 and S4. The files they replaced, `figures/figure3.png`,
`figures/figure_s2_fn_sensitivity.png` and `figures/figure_S4.png`, are kept so
that earlier versions stay reproducible.

## Four checks

| Check | Script |
|---|---|
| Share of training sites on label-homogeneous proteins, and separately on pure-positive proteins | `audit_ptm_benchmark.py`, `analysis/gate_homogeneity.py`, `analysis/benchmark_homogeneity.py` |
| Protein-identity baseline: score each site by its protein's training positive rate | `audit_ptm_benchmark.py` |
| Within-protein mean square ratio with and without the channel; null value 1.0 | `analysis/icc_by_task.py`, `analysis/icc_audit_v2.py` |
| Permutation control: does a permuted feature vector reproduce the gain? | `train_pdisjoint.py --cond shuffled` |

Recovery near 100% under permutation indicates a channel being used as a protein
identifier. A negative reading only shows the channel is *not* being used as an
identifier; it does not establish content use, since a pure noise channel
behaves the same way.

The evidence for transferable content is the paired comparison `real >
permuted`, and even that is conditional on how the test partition was built.
Under threshold-sampled negatives the real ESM-2 embedding beats the permutation
in 24 of 24 task-split pairs; under natural negatives it *loses* in six of eight
tasks and 17 of 24 pairs (corrected in v1.6, from `results/crosseval_verify.tsv`).
Report both test constructions.

Three of the four checks can be computed on an existing benchmark. The fourth —
the **natural-negative test**, evaluating the same trained model against
negatives drawn without the donor threshold — cannot. Under threshold-sampled
evaluation a frozen language model channel passes the permutation control in
every partition, shows no dose-response with homogeneity, and has a bootstrap
interval that excludes zero in all eight tasks on the positive side, while the
same trained models lose 0.13 AUPRC under natural negatives. Deciding whether a
protein-level channel generalises means rebuilding the negatives.

Confidence intervals throughout are protein-level
(`analysis/cluster_bootstrap.py` for the identity baseline,
`analysis/crosseval_bootstrap.py` for the cross-evaluation), not site-level; the
ICC values reported here contradict the independence assumption the latter makes.
For the cross-evaluation the current manuscript uses
`analysis/crosseval_bootstrap_v2.py`, which supersedes `crosseval_bootstrap.py`.

## Quick start

    python audit_ptm_benchmark.py --dir <dir with {ptm}_{train,test}.tsv>

Input files need the columns `protein`, `aa`, `pos`, `x` (sequence window), `y`.
Window length is not fixed; 31-mers and 61-mers both work.

To also report the identity baseline against a sequence model, supply the
model's AUROCs:

    python audit_ptm_benchmark.py --dir <dir> --seq-auroc seq_auroc.json

All scripts assume the repository root as the working directory. Paths default
to the environment in which the analysis was run and are overridable through
`PTM_AUDIT_BASE`; `rebuild_datasets.py` additionally reads `PTM_AUDIT_FASTA` and
`PTM_AUDIT_RETRAIN`, and `build_protein_features.py` and `node2vec_train.py`
take theirs as command-line arguments.

## Full pipeline

Dataset construction and training:

    rebuild_datasets.py       # rebuild datasets varying only the sampling rule
    make_pdisjoint_split.py   # CD-HIT clusters -> protein-disjoint partitions
    train_pdisjoint.py        # retrain under those partitions
    train_alldata.py          # same, on the authors' released protein-level split
    rerun_inference.py        # regenerate predictions from saved checkpoints

`train_pdisjoint.py` evaluates every trained ensemble on both reconstructions'
test partitions. The split file is shared, so `split == 'test'` selects the same
protein set either way; only which candidate sites appear as negatives differs.
Under a protein-disjoint split neither test set carries lookup leakage, since no
test protein appears in training. `--feat {ppi,esm,prott5}` selects the
protein-level channel.

Training runs, in order:

    sbatch slurm/train_pdisjoint_v2.slurm   # 144 cells, interaction embedding
    sbatch slurm/extract_esm.slurm          # language model feature extraction
    sbatch slurm/train_esm.slurm            # 96 cells, language model arm

Analysis:

    audit_ptm_benchmark.py    # dataset audit, T1-T6
    restricted_eval.py        # full vs. homogeneity-restricted subsets
    unseen_protein_eval.py    # evaluation on proteins absent from training
    bootstrap_ci.py           # paired bootstrap intervals (site-level; superseded)
    summarize_pdisjoint.py    # aggregate protein-disjoint runs
    make_table1.py            # assemble the main comparison
    make_figures.py           # Supplementary Figures S1 and S3

    analysis/gate_homogeneity.py     # homogeneity by construction and partition
    analysis/benchmark_homogeneity.py # the same statistics for any benchmark
    analysis/crosseval_summary.py    # 2x2x2 cross-evaluation
    analysis/crosseval_verify.py     # independent recomputation of the same cells
    analysis/crosseval_bootstrap.py  # paired cluster bootstrap, 192 cells
    analysis/merge_ci_shards.py      # collect the bootstrap array output
    analysis/depth_stratified.py     # within-task dose-response against depth
    analysis/mechanism_chain.py      # annotation-depth mechanism, end to end
    analysis/icc_by_task.py          # ICC and MSW, all tasks, both constructions
    analysis/icc_audit_v2.py         # the same, recomputed from the predictions
    analysis/cluster_bootstrap.py    # protein-level bootstrap
    analysis/homology_control.py     # nearest-neighbour separation, error-distance test
    analysis/sanity_posrate.py       # per-protein positive rate variability
    analysis/synthetic_reversal.py   # synthetic populations, demo and parameter sweep
    analysis/make_figures_v2.py      # Figures 2 and 4
    analysis/make_figure_merged.py   # Figure 3
    analysis/plot_depth_strata.py    # Supplementary Figure S4
    analysis/plot_perm_margin.py     # Supplementary Figure S5
    analysis/parse_musitedeep.py     # annotated FASTA -> (protein, pos, label)

`bootstrap_ci.py` resamples sites, which assumes site independence — contradicted
by the ICC values this work reports. Use `analysis/cluster_bootstrap.py`, which
resamples whole proteins and tests the paired difference of declines rather than
whether two marginal intervals overlap.

`analysis/patch_fig4.py` is a one-off that has already been applied;
`make_figure_merged.py` carries the fix. It is kept for the record because the
bug it corrects will bite anyone rerunning an earlier commit: feature families
were matched by `endswith('__esm')`, which was correct while only two families
existed and silently overwrote the interaction series once ProtT5 was added.

SLURM submission scripts are in `slurm/`.

## Reproducing the manuscript

| Manuscript item | Script |
|---|---|
| Table 1, homogeneity by construction | `audit_ptm_benchmark.py` section T3, cross-checked by `analysis/benchmark_homogeneity.py` |
| Tables 2 and 3, Supplementary Tables S3, S4, S11 | `analysis/crosseval_summary.py [--all]`, independently recomputed by `analysis/crosseval_verify.py` |
| Bootstrap intervals in Tables 2 and 3 | `slurm/crosseval_ci.slurm` then `analysis/merge_ci_shards.py` |
| Supplementary Tables S1, S2 | `analysis/cluster_bootstrap.py` |
| Supplementary Table S5, annotation depth | `analysis/section3_corrected.py` (`results/section3_20260920/ridge_summary.tsv`); `analysis/mechanism_chain.py` produced the version before the Section 3 revision |
| Supplementary Tables S6, S7, S8 | `audit_ptm_benchmark.py`, `restricted_eval.py` |
| Supplementary Table S9, current release | `train_alldata.py`, then `analysis/summarize_alldata.py` (AUROC, and the across-seed spread quoted in Note S6) and `analysis/s9_auprc_v2.py` (the AUPRC companion) |
| Supplementary Table S10, false-negative sensitivity | `analysis/fn_sensitivity.py`; the methylation cells rerun by `analysis/recompute_methyl_sensitivity.py` and merged by `analysis/assemble_section3_sensitivity.py`; then `analysis/make_table_s10.py --results results/section3_20260920` |
| Supplementary Table S12, cross-benchmark homogeneity | `analysis/benchmark_homogeneity.py` |
| Variance structure, all tasks | `analysis/icc_by_task.py [--feat esm]`, `analysis/icc_audit_v2.py` |
| Supplementary Note S2, homology control | `analysis/homology_control.py` |
| Supplementary Note S9, depth stratification | `analysis/depth_stratified.py` |
| Per-protein positive rate check | `analysis/sanity_posrate.py` |
| Share of deployment negatives below the threshold | `analysis/unobserved_weight.py` |

### The share of deployment negatives the threshold removes

`analysis/unobserved_weight.py` reports, per task, the fraction of the negatives
a predictor is actually applied to that sit on proteins carrying fewer than ten
annotated sites of the merged type — the region a threshold-sampled benchmark
supplies no negatives from at all. It weights each protein by the candidate
residues it contributes, so it needs only a residue count and an annotation
depth: no labels, no negatives from the benchmark, and no trained model. That
makes it computable for any released benchmark from its annotations plus a
reference proteome, and computable by a benchmark builder while choosing the
threshold rather than afterwards.

Two universes are reported because they answer different questions. `all` counts
every protein in the reference proteome, which is the population a predictor
meets — most proteins it sees carry no annotated site and therefore fall below
any positive threshold by construction. `seen` restricts to proteins already
carrying at least one site of the task, which is the population the
reconstructions cover and a conservative floor. The manuscript quotes both.

### Confidence intervals for the cross-evaluation

This describes the v1.5 procedure, which pooled the splits, so it is not the
statistic reported in Table 1. See `analysis/crosseval_bootstrap_v2.py` and
`docs/crosseval_bootstrap_v2.txt` for the procedure the current manuscript uses.

A paired cluster bootstrap over whole test proteins, 2,000 resamples for each of
the 192 cells. Within an iteration the baseline and the augmented model are
scored on an identical protein sample, so the interval is on their difference.
Each of the three partitions is resampled over its own proteins and the results
combined before scoring, so the interval carries partition variance as well.

    mkdir -p logs ci_shards
    sbatch slurm/crosseval_ci.slurm            # 192 array tasks, up to ~7 min each
    python analysis/merge_ci_shards.py ci_shards -o results/crosseval_ci.tsv

Seeds derive from the cell identity through sha256 rather than from position, so
a single shard reproduces exactly the row a whole-set run would produce, and a
failed task can be resubmitted alone with `sbatch --array=N slurm/crosseval_ci.slurm`.
Without Slurm, `python analysis/crosseval_bootstrap.py --root <runs>` does the
whole set serially in a few hours.

### Within-task dose-response against annotation depth

Bins the test proteins of each task by annotation depth, measured on the
unrestricted reconstruction, and reports the change in each bin's mean position
in the score ranking of the whole task. Depth bins are crossed with
candidate-site-count bins because depth and protein size are correlated.

    python analysis/depth_stratified.py --root <runs> --data rebuilt \
        --family esm2 --train replica --eval rebuilt -o results/depth_esm2_thr.tsv
    python analysis/depth_stratified.py --root <runs> --data rebuilt \
        --family esm2 --train rebuilt --eval rebuilt -o results/depth_esm2_unr.tsv
    python analysis/plot_depth_strata.py results/depth_esm2_thr.tsv \
        results/depth_esm2_unr.tsv -o figures/figure_S4.png

The `--train rebuilt` run is the negative control: the gradient should not
survive it. Repeat with `--family interaction` and `--family prott5`.

The primary quantity is the mean score percentile, not a within-bin AUPRC. The
shortcut acts between proteins — a model that has learned that shallow
annotation implies a positive label raises every site of such a protein by a
similar amount — so inside a bin of uniformly shallow proteins the ranking
barely moves, and a within-bin metric cannot see the effect it is meant to
measure. Within-bin AUPRC is in the output table and is close to noise for that
reason.

### Label homogeneity across benchmarks

    python analysis/benchmark_homogeneity.py --data retrain \
        --name "DeepMVP earlier release" --only train -o results/cross_benchmark.tsv
    python analysis/benchmark_homogeneity.py --data rebuilt \
        --name "This work, unrestricted" --only all -o results/cross_benchmark.tsv

For MusiteDeep, whose releases are full-length sequences with modified residues
marked by `#`, convert first:

    python analysis/parse_musitedeep.py \
        --root <MusiteDeep_web>/MusiteDeep/testdata --only train \
        -o musitedeep_sites.tsv
    python analysis/benchmark_homogeneity.py --data musitedeep_sites.tsv \
        --name MusiteDeep -o results/cross_benchmark.tsv

Run it on `retrain` before running it on anything else. It must reproduce the T3
output of `audit_ptm_benchmark.py`; if it does not, the columns have been read
wrongly and the figures for other benchmarks mean nothing. That check is how the
homogeneous-site column of Table 1 was found to be wrong for two tasks.

### Independent recomputation

`analysis/crosseval_verify.py` recomputes every cell from the per-site
predictions without calling `crosseval_summary.py`. Rerunning the original
script would only show that it is deterministic; a second implementation can
disagree. The two agree to four decimal places on all 192 cells.

`analysis/icc_audit_v2.py` stands in the same relation to `icc_by_task.py`.
Ratios there are formed after averaging over partitions, never by averaging
per-partition ratios: ICC(1,1) can be near zero, and a mean of per-partition
folds then diverges.

### Reproduction control

`pdisjoint_runs/` holds an earlier round and is retained as a reproduction
control: across the 120 cells shared with the current round the mean difference
in AUROC is −0.0000 (sd 0.0033, maximum 0.0171, none above 0.02, uniform across
tasks). Exact reproduction is not possible because cuDNN provides no
deterministic GRU implementation, so the manuscript reports the second round
throughout and the first is kept for comparison.

Predictions recorded before each experiment was run, including the two that
turned out wrong, are in `docs/docs_prediction_crosseval.txt`.

### A note on the merged methylation task

The file released as methylation K contains both lysine and arginine sites. Our
reconstructions rebuild each task from the same positive annotations, so the
overlap with the methylation R task is inherited rather than removed, and it is
labelled `Methylation K/R` in every table. Counts over the eight tasks are
therefore not counts over eight independent units; collapsing the two shifts the
reported means by 4% to 17% without changing any sign.

## Embeddings

    node2vec_train.py         # STRING v12.0 physical subnetwork -> ENSP embeddings
    (UniProt ID mapping web service -> idmapping_*.tsv)
    map_string_to_uniprot.py  # emb + idmapping -> node2vec_with_uniprot.csv
    build_protein_features.py # csv -> protein_features_ppi.npy + protein_ids_ppi.json
    extract_esm.py            # frozen mean-pooled ESM-2 650M -> protein_features_esm.npy
    extract_prott5.py         # frozen mean-pooled ProtT5-XL-U50

The ENSP-to-UniProt lookup itself goes through the UniProt web mapping service
and is not scriptable; its output table is archived in `data/` so the chain can
be verified end to end. Collapsing that table onto the embedding *is* scripted,
in `map_string_to_uniprot.py`: 377 of the 18,423 STRING nodes carry no
accession and 37 accessions are reached by more than one node, so a
duplicate-handling policy has to be chosen. Run it with `--verify` first — it
rebuilds the matrix under each policy and reports which one reproduces the
archived `protein_features_ppi.npy` — then rerun with that policy. This is the
step whose arithmetic Supplementary Note S10 reports. Because node2vec is
stochastic, rerunning
`node2vec_train.py` produces a functionally equivalent but not bitwise-identical
embedding — the archived matrix is the one used in the manuscript.

`(p, q)` was selected from {(1, 0.5), (1, 1), (4, 1), (4, 4)} by maximising the
mean cosine similarity of connected node pairs, which chose (1, 0.5). That
criterion would also reward a degenerate embedding in which all vectors
coincide, so we checked: the mean cosine similarity of randomly paired proteins
is 0.186, and no collapse occurred.

That check does not transfer to a protein language model. Frozen mean-pooled
ESM-2 embeddings have a raw mean pairwise cosine near 0.87, which reflects the
known anisotropy of transformer representations rather than degeneracy: after
subtracting the mean vector it falls to approximately zero, and the effective
rank is 200.5 of 1280. Report the centred cosine and the effective rank for pLM
features, not the raw cosine. The features are used uncentred, so that the two
channels differ only in their source.

## Synthetic populations

`analysis/synthetic_reversal.py` generates data satisfying only the three
structural conditions — items grouped into entities, negatives drawn only from
entities carrying at least T positives, and an entity-constant feature encoding
that count — with no reference to sequence, structure or biology, and fits a
plain multilayer perceptron rather than the architecture used on the real data.

The single-point comparison and the parameter scan use different replicate
counts, and neither is the script default. The manuscript reports:

    python analysis/synthetic_reversal.py --mode demo  --seeds 5 --out results
    python analysis/synthetic_reversal.py --mode sweep --seeds 3 --out results

`--seeds` defaults to 3. Running `demo` without `--seeds 5` produces a different
table and overwrites `results/synthetic_demo.tsv`, which is how that file was
lost once; the version in the repository is the five-seed run the manuscript
quotes.

The sweep locates the sign change between T = 5 and T = 10, shows that
thresholds of one and two are indistinguishable from no threshold at all, and
shows that under a permissive rule the contribution grows more positive as the
feature encodes depth better while above the boundary the same gradient runs the
other way. The fidelity axis reports the measured cross-validated ridge R² from
the feature to log depth, not the `r2_target` parameter, which is roughly twice
as large.

## Data

`data/` holds the interaction embeddings, the STRING-to-UniProt mapping table
they derive from, and the protein-disjoint partition assignments. The two
analysed PTM releases, their download dates and their checksums, and the
provenance of the language model embeddings, are described in `DATA.md`; the
releases themselves are archived separately because the individually posted
files on the source portal are expected to change.

Reconstructed datasets, partition assignments, per-site predictions for every
cell of the cross-evaluation, trained model weights and all three feature
matrices are deposited in a Zenodo data record, concept DOI
10.5281/zenodo.21670042, which resolves to its latest version. That version also
holds the ratio-matched controls and the mouse dataset. This repository is
archived under concept DOI 10.5281/zenodo.21993819.

### External benchmarks

MusiteDeep, from `github.com/duolinwang/MusiteDeep_web`, training files under
`MusiteDeep/testdata/<modification>/<residues>/`. Full-length sequences with
modified residues marked by `#`; every unmarked residue of the target type is a
negative, so no donor threshold applies.

PhosIDN, from `github.com/ustchangyuanyang/PhosIDN`. Code and weights only.
`methods/dataprocess_train.py` reads labels from column 0 of a CSV the
repository does not contain, so its negative construction cannot be checked from
what is public.

Neither clone is included here; both are large and belong to their authors.

## Requirements

Analysis and training:

    python >= 3.9, numpy, pandas, scipy, scikit-learn, matplotlib
    torch (training only)
    CD-HIT >= 4.8 (protein-disjoint partitioning)

Embedding generation only:

    networkx, node2vec, gensim   # interaction embeddings
    fair-esm, transformers       # ESM-2 and ProtT5 embeddings

Pinned versions are in `requirements.txt`.

## Citation

[to be completed]

## Licence

MIT. See `LICENSE`.

### SLURM scripts

    slurm/cdhit.slurm                # CD-HIT clustering for the protein-disjoint split
    slurm/train_pdisjoint_v2.slurm   # 144 cells, interaction embedding, both test sets
    slurm/extract_esm.slurm          # language model feature extraction
    slurm/train_esm.slurm            # 96 cells, language model arm
    slurm/alldata.slurm              # training on the authors' current release
    slurm/rerun_inference.slurm      # regenerate published-protocol predictions
    slurm/crosseval_ci.slurm         # 192-cell paired cluster bootstrap

`train_pdisjoint.slurm` and `shuffle_pdisjoint.slurm` are the earlier round,
retained because `pdisjoint_runs/` was produced by them and is kept as a
reproduction control. `shuffle.slurm` is the permutation control under the
published protocol.

### Archived outputs

`results/` holds the console output of every analysis script, so the numbers in
the manuscript can be checked without rerunning anything:

    crosseval_ppi.txt / crosseval_esm.txt   the 2x2x2 cross-evaluation
    crosseval_verify.tsv                    the same, independently recomputed
    crosseval_ci.tsv                        bootstrap intervals, 192 cells (pooled; superseded)
    crosseval_ci_v2.tsv                     bootstrap intervals, the statistic in Table 1
    table1_cells.tsv                        Table 1 of the current manuscript
    mechanism_ppi.txt / mechanism_esm.txt   the annotation-depth chain
    icc_ppi.txt / icc_esm.txt               ICC and within-protein mean squares
    icc_audit.tsv                           the same, recomputed from predictions
    cluster_bootstrap_out.txt               protein-level intervals
    homology_control.txt                    nearest-neighbour separation
    sanity_posrate.txt                      per-protein positive rate spread
    depth_*_thr.tsv / depth_*_unr.tsv       depth stratification, three channels
    cross_benchmark.tsv                     label homogeneity across benchmarks
    synthetic_demo.tsv                      the five-seed single-point comparison
    figures/figure6_sweep.tsv               the full synthetic parameter sweep

`table1.txt`, `summary.json`, `unseen_eval.txt` and the two checksum files are
from the earlier round.

## Section 3 revision (2026-09-20)

The revised Section 3 and the revised false-negative sensitivity analysis
were computed with three scripts that read the existing reconstructions and
prediction files only; no networks are retrained.

    analysis/section3_corrected.py             annotation-depth probe with size
                                               and length controls and residual R2,
                                               label-rate and score-depth
                                               correlations, fixed-neighbour
                                               analysis, depth by candidate-count
                                               strata; writes Figure 2 and
                                               Supplementary Figure S4
    analysis/recompute_methyl_sensitivity.py   reruns the false-negative
                                               sensitivity cells affected by
                                               de-duplicating methylation positions,
                                               one shard per run
    analysis/assemble_section3_sensitivity.py  merges those cells with the
                                               unchanged tasks into
                                               fn_sensitivity_{ppi,esm}.tsv

Outputs are in `results/section3_20260920/`. `manifest.json` records the
SHA-256 of the `section3_corrected.py` that produced them and the definitions
used (depth: unique protein-position annotations per merged type; size: full
target-residue count from the reference FASTA).

Paths. The scripts are committed as they were run on the UConn HPC, where the
working root is `/home/FCAM/juli/HRP` (it holds this repository, the
reconstructions and the prediction files; `manifest.json` records it as
`input_root`). `section3_corrected.py` and `recompute_methyl_sensitivity.py`
take that root as `--root` and an output directory as `--out`;
`recompute_methyl_sensitivity.py` also takes `--shard`, and shards 0 to 11 were
run, writing to `section3_revision_20260920/fn_shards/shard00` to `shard11`.
`assemble_section3_sensitivity.py` has the root written into it: it reads
`/home/FCAM/juli/HRP/ptm-audit/results/fn_sensitivity_{ppi,esm}.tsv` and the
shard outputs under `/home/FCAM/juli/HRP/section3_revision_20260920/`, and
writes to `/home/FCAM/juli/HRP/section3_revision_20260920/results/`. Change
`root` at the top of that script to run it elsewhere.

Figure 2 and Supplementary Figure S4 as they appear in the manuscript were
redrawn by `analysis/plot_corrected.py` from six of the tables that
`section3_corrected.py` writes (`ridge_acetylation_oof.tsv`,
`ridge_summary.tsv`, `score_depth_summary.tsv`, `knn_summary.tsv`,
`depth_plot_splits.tsv`, `depth_by_size_splits.tsv`). The two PNGs in
`results/section3_20260920/` are those redrawn figures; the earlier renderings
written by `section3_corrected.py` itself are in commit `13b1817`, and rerunning
that script writes the earlier layout again.

`plot_corrected.py` is committed as it was run, from a separate working folder
on a laptop. It reads its tables from, and writes the two PNGs to, the folder
`results/` one level above the script (`P` near the top of the script); inside
this repository that would be the top-level `results/`, so set `P` to
`results/section3_20260920` to run it here. It also puts a folder `render_deps`,
one level above the script, at the front of `sys.path`; that folder held a local
copy of matplotlib and is not needed where matplotlib is installed. With the
committed tables and matplotlib 3.11, the script reproduces both PNGs pixel for
pixel.
