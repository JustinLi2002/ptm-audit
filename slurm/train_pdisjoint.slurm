#!/bin/bash
### One array task per (ptm, dataset, condition, split seed) = 96 cells.
### Each task trains --n-models models and evaluates on the common test set.
#SBATCH --job-name=pdisj
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=108G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jli89@tulane.edu
#SBATCH --array=0-95%8
#SBATCH -e logs/pdisj_%A_%a.err
#SBATCH -o logs/pdisj_%A_%a.log

source /home/FCAM/juli/HRP/my-venv/bin/activate
cd /home/FCAM/juli/HRP
mkdir -p logs pdisjoint_runs

PTMS=(acetylation_k glycosylation_n methylation_k methylation_r \
      phosphorylation_st phosphorylation_y sumoylation_k ubiquitination_k)
DATASETS=(replica rebuilt)
CONDS=(baseline ppi)
SPLITS=(0 1 2)

i=$SLURM_ARRAY_TASK_ID
s=$(( i % 3 ));            i=$(( i / 3 ))
c=$(( i % 2 ));            i=$(( i / 2 ))
d=$(( i % 2 ));            i=$(( i / 2 ))
p=$(( i % 8 ))

PTM=${PTMS[$p]}; DS=${DATASETS[$d]}; CD=${CONDS[$c]}; SP=${SPLITS[$s]}
OUT="pdisjoint_runs/${PTM}__${DS}__${CD}__split${SP}.json"

if [ -f "$OUT" ]; then
    echo "already done: $OUT"; exit 0
fi

echo "task $SLURM_ARRAY_TASK_ID -> $PTM / $DS / $CD / split $SP"
echo "Host: $(hostname)  Time: $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

python train_pdisjoint.py --ptm "$PTM" --dataset "$DS" --cond "$CD" \
       --split-seed "$SP" --n-models 2

echo "Done at $(date)"
