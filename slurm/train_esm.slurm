#!/bin/bash
### E3: frozen mean-pooled ESM-2 650M as the protein-level channel.
### 96 cells = 8 ptm x 2 dataset x 2 cond x 3 split.
### Baselines are reused from the PPI runs -- they contain no vector at all.
### A frozen pLM embedding is constant within a protein, structurally identical
### to the interaction embedding, and encodes annotation depth 2-3x better
### (CV R2 .277-.575 vs .126-.235), so the same 2x2 should appear more strongly.
#SBATCH --job-name=esm3
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jli89@tulane.edu
#SBATCH --array=0-95%8
#SBATCH -e logs/esm3_%A_%a.err
#SBATCH -o logs/esm3_%A_%a.log

source /home/FCAM/juli/HRP/my-venv/bin/activate
cd /home/FCAM/juli/HRP
mkdir -p logs pdisjoint_runs_v2 checkpoints_pdisjoint

PTMS=(acetylation_k glycosylation_n methylation_k methylation_r \
      phosphorylation_st phosphorylation_y sumoylation_k ubiquitination_k)
DATASETS=(replica rebuilt)
CONDS=(ppi shuffled)
SPLITS=(0 1 2)

i=$SLURM_ARRAY_TASK_ID
s=$(( i % 3 )); i=$(( i / 3 ))
c=$(( i % 2 )); i=$(( i / 2 ))
d=$(( i % 2 )); i=$(( i / 2 ))
p=$(( i % 8 ))

PTM=${PTMS[$p]}; DS=${DATASETS[$d]}; CD=${CONDS[$c]}; SP=${SPLITS[$s]}
OUT="pdisjoint_runs_v2/${PTM}__${DS}__${CD}__split${SP}__esm.json"
[ -f "$OUT" ] && { echo "already done: $OUT"; exit 0; }

echo "task $SLURM_ARRAY_TASK_ID -> $PTM / $DS / $CD / split $SP / esm"
echo "Host: $(hostname)  Time: $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

python ptm-audit/train_pdisjoint.py --ptm "$PTM" --dataset "$DS" --cond "$CD" \
       --split-seed "$SP" --feat esm --n-models 2

echo "Done at $(date)"
