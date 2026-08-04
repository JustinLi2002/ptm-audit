#!/bin/bash
### 144 cells = 8 ptm x 2 dataset x 3 cond x 3 split.
### 24 of these (replica x shuffled) were never run before.
### Each cell evaluates on BOTH test constructions in one pass.
#SBATCH --job-name=pdisj2
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=68G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jli89@tulane.edu
#SBATCH --array=0-143%8
#SBATCH -e logs/pdisj2_%A_%a.err
#SBATCH -o logs/pdisj2_%A_%a.log

source /home/FCAM/juli/HRP/my-venv/bin/activate
cd /home/FCAM/juli/HRP
mkdir -p logs pdisjoint_runs_v2 checkpoints_pdisjoint

PTMS=(acetylation_k glycosylation_n methylation_k methylation_r \
      phosphorylation_st phosphorylation_y sumoylation_k ubiquitination_k)
DATASETS=(replica rebuilt)
CONDS=(baseline ppi shuffled)
SPLITS=(0 1 2)

i=$SLURM_ARRAY_TASK_ID
s=$(( i % 3 )); i=$(( i / 3 ))
c=$(( i % 3 )); i=$(( i / 3 ))
d=$(( i % 2 )); i=$(( i / 2 ))
p=$(( i % 8 ))

PTM=${PTMS[$p]}; DS=${DATASETS[$d]}; CD=${CONDS[$c]}; SP=${SPLITS[$s]}
OUT="pdisjoint_runs_v2/${PTM}__${DS}__${CD}__split${SP}.json"
[ -f "$OUT" ] && { echo "already done: $OUT"; exit 0; }

echo "task $SLURM_ARRAY_TASK_ID -> $PTM / $DS / $CD / split $SP"
echo "Host: $(hostname)  Time: $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

python ptm-audit/train_pdisjoint.py --ptm "$PTM" --dataset "$DS" --cond "$CD" \
       --split-seed "$SP" --n-models 2

echo "Done at $(date)"
