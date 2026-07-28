# ptm-audit

Diagnostics for label leakage through protein-level feature channels in PTM
site prediction. Companion code for [manuscript in preparation].

## What the leakage is

Negative samples in PTM benchmarks are drawn only from proteins above a
site-count threshold. Proteins below it contribute positives and no negatives,
so every one of their sites carries the same label. Any feature that is constant
within a protein — an interaction embedding, a kinase-substrate prior — can then
act as a lookup key, provided the same protein appears in both partitions.

## Three checks

| Script | Question |
|---|---|
| `audit_ptm_benchmark.py` | How much of the dataset is label-homogeneous, and how much of the test set does that determine? |
| `bootstrap_ci.py` | Does a protein-identity-only classifier approach the model? |
| `train_pdisjoint.py --cond shuffled` | Does a permuted feature vector reproduce the gain? |

Recovery near 100% under permutation indicates a channel used as a protein
identifier. Recovery at or below zero indicates a channel used for its content.

## Quick start

    python audit_ptm_benchmark.py --dir <dir with {ptm}_{train,test}.tsv>

Input files need columns `protein`, `aa`, `pos`, `x` (sequence window), `y`.

## Full pipeline

    rebuild_datasets.py      # rebuild datasets varying only the sampling rule
    make_pdisjoint_split.py  # CD-HIT clusters -> protein-disjoint partitions
    train_pdisjoint.py       # retrain under those partitions
    train_alldata.py         # same, on the authors' released protein-level split
    rerun_inference.py       # regenerate predictions from saved checkpoints
    make_table1.py           # assemble the main comparison
    bootstrap_ci.py          # paired bootstrap intervals

## Requirements

python >= 3.9, numpy, pandas, scikit-learn, torch (training only), CD-HIT >= 4.8

## Data

Analysed releases and download dates are listed in `DATA.md`.
