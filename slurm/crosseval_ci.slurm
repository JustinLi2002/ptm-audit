#!/bin/bash
#SBATCH --job-name=xeval_ci
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A100:1
#SBATCH --array=0-191%32
#SBATCH --mem=20G
#SBATCH --time=08:00:00
#SBATCH --output=logs/xeval_ci_%A_%a.out
#SBATCH --error=logs/xeval_ci_%A_%a.err
#
# One array task per (family, training construction, evaluation construction,
# task) cell: 6 x 2 x 2 x 8 = 192.
#
#   mkdir -p logs ci_shards
#   sbatch crosseval_ci.slurm
#
# Each shard writes its own row, so a failed task can be resubmitted alone:
#   sbatch --array=57 crosseval_ci.slurm
#
# Seeds derive from the cell identity through sha256, not from position, so a
# shard reproduces exactly the row a whole-set run would have produced and
# reruns are identical.

set -euo pipefail

ROOT="${ROOT:-$HOME/HRP/pdisjoint_runs_v2}"
OUTDIR="${OUTDIR:-$HOME/HRP/ci_shards}"
B="${B:-2000}"

mkdir -p "$OUTDIR"

# BLAS threads fight with the array's own parallelism; one thread per task is
# faster here because the work is many small metric evaluations, not big matmuls.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

python "$HOME/HRP/crosseval_bootstrap.py" \
    --root "$ROOT" \
    -B "$B" \
    --shard "$SLURM_ARRAY_TASK_ID" \
    -o "$OUTDIR/cell_$(printf '%03d' "$SLURM_ARRAY_TASK_ID").tsv"
