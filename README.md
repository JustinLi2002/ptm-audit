# ptm-audit

Diagnostics for label leakage through protein-level feature channels in PTM site
prediction. Companion code for [manuscript in preparation].

## What the leakage is

Negative samples in PTM benchmarks are drawn only from proteins carrying at
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

## Four checks

| Check | Script |
|---|---|
| Share of training sites on label-homogeneous proteins, and separately on pure-positive proteins | `audit_ptm_benchmark.py`, `analysis/gate_homogeneity.py` |
| Protein-identity baseline: score each site by its protein's training positive rate | `audit_ptm_benchmark.py` |
| Within-protein mean square ratio with and without the channel; null value 1.0 | `analysis/icc_by_task.py` |
| Permutation control: does a permuted feature vector reproduce the gain? | `train_pdisjoint.py --cond shuffled` |

Recovery near 100% under permutation indicates a channel being used as a protein
identifier. A negative reading only shows the channel is *not* being used as an
identifier; it does not establish content use, since a pure noise channel
behaves the same way.

The evidence for transferable content is the paired comparison `real >
permuted`, and even that is conditional on how the test partition was built.
Under threshold-sampled negatives the real embedding beats the permutation in
8/8 tasks; under natural negatives it *loses* in the five tasks where the
shortcut is available. Report both test constructions.

Three of the four checks can be computed on an existing benchmark. The fourth
cannot: under threshold-sampled evaluation a frozen language model channel
passes the permutation control in every partition, shows no dose-response with
homogeneity, and improves AUPRC in every task, while the same trained models
lose 0.13 AUPRC under natural negatives. Deciding whether a protein-level
channel generalises requires a test partition whose negatives were not drawn
under the rule being tested, which in practice means rebuilding them.

Confidence intervals throughout are protein-level (`analysis/cluster_bootstrap.py`),
not site-level; the ICC values reported here contradict the independence
assumption the latter makes.

## Quick start

    python audit_ptm_benchmark.py --dir <dir with {ptm}_{train,test}.tsv>

Input files need the columns `protein`, `aa`, `pos`, `x` (sequence window), `y`.
Window length is not fixed; 31-mers and 61-mers both work.

To also report the identity baseline against a sequence model, supply the
model's AUROCs:

    python audit_ptm_benchmark.py --dir <dir> --seq-auroc seq_auroc.json

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
test protein appears in training. `--feat {ppi,esm}` selects the protein-level
channel.

Analysis:

    audit_ptm_benchmark.py    # dataset audit, T1-T6
    restricted_eval.py        # full vs. homogeneity-restricted subsets
    unseen_protein_eval.py    # evaluation on proteins absent from training
    bootstrap_ci.py           # paired bootstrap intervals (site-level; superseded)
    summarize_pdisjoint.py    # aggregate protein-disjoint runs
    make_table1.py            # assemble the main comparison
    make_figures.py           # Figures 2 and 3

    analysis/gate_homogeneity.py    # homogeneity by construction and partition
    analysis/crosseval_summary.py   # 2x2x2 cross-evaluation
    analysis/mechanism_chain.py     # annotation-depth mechanism, end to end
    analysis/icc_by_task.py         # ICC and MSW, all tasks, both constructions
    analysis/cluster_bootstrap.py   # protein-level bootstrap for T2/S1/S2
    analysis/homology_control.py    # nearest-neighbour separation, error-distance test
    analysis/sanity_posrate.py      # per-protein positive rate variability
    analysis/make_figures_v2.py     # Figures 4-7, from the run outputs

`bootstrap_ci.py` resamples sites, which assumes site independence — contradicted
by the ICC values this work reports. Use `analysis/cluster_bootstrap.py`, which
resamples whole proteins and tests the paired difference of declines rather than
whether two marginal intervals overlap.

SLURM submission scripts for the cluster runs are in `slurm/`.

Figure 1 is a hand-drawn schematic (`figures/fig1_mechanism.svg`). Figures 2 and
3 are produced by `make_figures.py` from hard-coded values. Figures 4 to 7 are
produced by `analysis/make_figures_v2.py`, which reads `pdisjoint_runs_v2/` and
the feature matrices directly, so every panel is reproducible from the archived
data.

## Reproducing the manuscript

| Manuscript item | Script |
|---|---|
| Table 1, homogeneity by construction | `analysis/gate_homogeneity.py` |
| Tables 2, S1, S2 | `analysis/cluster_bootstrap.py` |
| Tables 3, 4 | `audit_ptm_benchmark.py`, `restricted_eval.py` |
| Tables 5, 6, S3, S4 | `analysis/crosseval_summary.py [--feat esm]` |
| Table 7, current release | `train_alldata.py`, `summarize_pdisjoint.py` |
| Table S5, annotation depth | `analysis/mechanism_chain.py` |
| Figures 2, 3 | `make_figures.py` |
| Figures 4-7 | `analysis/make_figures_v2.py` |
| Variance structure, all tasks | `analysis/icc_by_task.py [--feat esm]` |
| Homology proximity control | `analysis/homology_control.py` |
| Per-protein positive rate check | `analysis/sanity_posrate.py` |

Training runs, in order:

    sbatch slurm/train_pdisjoint_v2.slurm   # 144 cells, interaction embedding
    sbatch slurm/extract_esm.slurm          # language model feature extraction
    sbatch slurm/train_esm.slurm            # 96 cells, language model arm

`pdisjoint_runs/` holds an earlier round and is retained as a reproduction
control: across the 120 shared cells the mean difference in AUROC between the
two rounds is −0.0000 (sd 0.0033, max 0.0171, none above 0.02, uniform across
tasks). Exact reproduction is not possible because cuDNN provides no
deterministic GRU implementation, so the manuscript reports the second round
throughout and the first is kept for comparison.

Predictions recorded before each experiment was run, including the two that
turned out wrong, are in `docs/docs_prediction_crosseval.txt`.

## Embeddings

    node2vec_train.py         # STRING v12.0 physical subnetwork -> ENSP embeddings
    (UniProt ID mapping web service -> node2vec_with_uniprot.csv)
    build_protein_features.py # csv -> protein_features_ppi.npy + protein_ids_ppi.json
    extract_esm.py            # frozen mean-pooled ESM-2 650M -> protein_features_esm.npy

The ENSP-to-UniProt step was performed through the UniProt web mapping service
and is not scriptable; the resulting table is archived in `data/` so the chain
can be verified end to end. Because node2vec is stochastic, rerunning
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

## Data

`data/` holds the interaction embeddings, the STRING-to-UniProt mapping table
they derive from, and the protein-disjoint partition assignments. The two
analysed PTM releases, their download dates and their checksums, and the
provenance of the language model embeddings, are described in `DATA.md`; the
releases themselves are archived separately because the individually posted
files on the source portal are expected to change.

Reconstructed datasets, partition assignments, per-site predictions for every
cell of the cross-evaluation, and both feature matrices are deposited at
[Zenodo DOI to be added].

## Requirements

Analysis and training:

    python >= 3.9, numpy, pandas, scipy, scikit-learn, matplotlib
    torch (training only)
    CD-HIT >= 4.8 (protein-disjoint partitioning)

Embedding generation only:

    networkx, node2vec, gensim   # interaction embeddings
    fair-esm                     # ESM-2 embeddings

Pinned versions are in `requirements.txt`.

## Citation

[to be completed]

## Licence

MIT. See `LICENSE`.
