# ptm-audit

Diagnostics for label leakage through protein-level feature channels in PTM site
prediction. Companion code for [manuscript in preparation].

## What the leakage is

Negative samples in PTM benchmarks are drawn only from proteins carrying at
least a threshold number of annotated sites. Proteins below the threshold
contribute positives and no negatives, so every one of their sites carries the
same label. Any feature that is constant within a protein — an interaction
embedding, a kinase-substrate prior — can then act as a key into a lookup table,
provided the same protein appears in both the training and test partitions.

Two conditions must therefore hold for the channel to be exploitable: a
permissive sampling rule, and a split that places the same protein on both
sides. Correcting either one closes it. Correcting only the split still leaves
the training set distorted.

## Three checks

| Script | Question |
|---|---|
| `audit_ptm_benchmark.py` | How much of the dataset is label-homogeneous, how much of the test set does that determine, and how far does a protein-identity-only classifier get? (tables T3, T4, T5) |
| `train_pdisjoint.py --cond shuffled` | Does a permuted feature vector reproduce the gain? |
| `bootstrap_ci.py` | Paired bootstrap intervals for the two preceding. |

Recovery near 100% under permutation indicates a channel being used as a protein
identifier. Recovery at or below zero indicates a channel being used for its
content.

## Quick start

    python audit_ptm_benchmark.py --dir <dir with {ptm}_{train,test}.tsv>

Input files need the columns `protein`, `aa`, `pos`, `x` (sequence window), `y`.
Window length is not fixed; 31-mers and 61-mers both work.

To also report the identity baseline against a sequence model, supply the
model's AUROCs:

    python audit_ptm_benchmark.py --dir <dir> --seq-auroc seq_auroc.json

## Full pipeline

Dataset construction and evaluation:

    rebuild_datasets.py       # rebuild datasets varying only the sampling rule
    make_pdisjoint_split.py   # CD-HIT clusters -> protein-disjoint partitions
    train_pdisjoint.py        # retrain under those partitions
    train_alldata.py          # same, on the authors' released protein-level split
    rerun_inference.py        # regenerate predictions from saved checkpoints

Analysis:

    audit_ptm_benchmark.py    # dataset audit, T1-T6
    restricted_eval.py        # full vs. homogeneity-restricted subsets
    unseen_protein_eval.py    # evaluation on proteins absent from training
    bootstrap_ci.py           # paired bootstrap intervals
    summarize_pdisjoint.py    # aggregate protein-disjoint runs
    make_table1.py            # assemble the main comparison
    make_figures.py           # Figures 2, 3 and 5

SLURM submission scripts for the cluster runs are in `slurm/`.

Figures 2, 3 and 5 are produced by `make_figures.py` from hard-coded values.
Figure 1 is a hand-drawn schematic (`figures/fig1_mechanism.svg`). Figure 4 is
generated from per-site predictions.

## Embeddings

    node2vec_train.py         # STRING v12.0 physical subnetwork -> ENSP embeddings
    (UniProt ID mapping web service -> node2vec_with_uniprot.csv)
    build_protein_features.py # csv -> protein_features_ppi.npy + protein_ids_ppi.json

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

## Data

`data/` holds the interaction embeddings, the STRING-to-UniProt mapping table
they derive from, and the protein-disjoint partition assignments. The two
analysed PTM releases, their download dates and their checksums are described in
`DATA.md`; the releases themselves are archived separately because the
individually posted files on the source portal are expected to change.

## Requirements

Analysis and training:

    python >= 3.9, numpy, pandas, scikit-learn, matplotlib
    torch (training only)
    CD-HIT >= 4.8 (protein-disjoint partitioning)

Embedding generation only:

    networkx, node2vec, gensim

## Citation

[to be completed]

## Licence

MIT. See `LICENSE`.
