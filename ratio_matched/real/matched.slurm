#!/bin/bash
### v4 ratio-matched control (PREREG_real.md).
### 18 cells = 2 ptm x 3 cond x 3 split, matched training construction, ESM-2; two inits each.
### Each cell evaluates on BOTH existing test constructions in one pass.
#SBATCH --job-name=matched4
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jli89@tulane.edu
#SBATCH --array=0-17%8
#SBATCH -e logs/matched4_%A_%a.err
#SBATCH -o logs/matched4_%A_%a.log

source /home/FCAM/juli/HRP/my-venv/bin/activate
cd /home/FCAM/juli/HRP
mkdir -p logs pdisjoint_runs_matched_v4 checkpoints_matched_v4

PTMS=(acetylation_k phosphorylation_st)
CONDS=(baseline ppi shuffled)
SPLITS=(0 1 2)

i=$SLURM_ARRAY_TASK_ID
s=$(( i % 3 )); i=$(( i / 3 ))
c=$(( i % 3 )); i=$(( i / 3 ))
p=$(( i % 2 ))

PTM=${PTMS[$p]}; CD=${CONDS[$c]}; SP=${SPLITS[$s]}
FEAT=esm; SUFFIX=__esm
if [ "$CD" = baseline ]; then FEAT=ppi; SUFFIX=; fi
OUT="pdisjoint_runs_matched_v4/${PTM}__matched__${CD}__split${SP}${SUFFIX}.json"
[ -f "$OUT" ] && { echo "already done: $OUT"; exit 0; }

echo "task $SLURM_ARRAY_TASK_ID -> $PTM / matched / $CD / split $SP / $FEAT"
echo "Host: $(hostname)  Time: $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

python ratio_matched_v4/train_matched.py --ptm "$PTM" --dataset matched --cond "$CD" \
       --split-seed "$SP" --feat "$FEAT" --n-models 2

echo "Done at $(date)"
