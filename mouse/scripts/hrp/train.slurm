#!/bin/bash
### Mouse replication (PREREG_mouse.md): 27 cells = 3 training constructions x 3 conditions x 3 splits, ESM-2,
### two initialisations each (54 fits). Each cell is scored on both test constructions in one pass. Cells with a
### result are skipped. The GPU is checked before training, and a failed run fails the job.
#SBATCH --job-name=mouse_train
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --exclude=mantis-048,mantis-075
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jli89@tulane.edu
#SBATCH --array=0-26%9
#SBATCH -e /home/FCAM/juli/HRP/mouse_v1/logs/train_%A_%a.err
#SBATCH -o /home/FCAM/juli/HRP/mouse_v1/logs/train_%A_%a.log
set -o pipefail
export PTM_AUDIT_BASE=/home/FCAM/juli/HRP/mouse_v1
source /home/FCAM/juli/HRP/my-venv/bin/activate
cd "$PTM_AUDIT_BASE"

DATASETS=(replica rebuilt matched)
CONDS=(baseline ppi shuffled)
i=$SLURM_ARRAY_TASK_ID
s=$(( i % 3 )); i=$(( i / 3 ))
c=$(( i % 3 )); i=$(( i / 3 ))
d=$(( i % 3 ))
DS=${DATASETS[$d]}; CD=${CONDS[$c]}; SP=$s
FEAT=esm; SUFFIX=__esm
if [ "$CD" = baseline ]; then FEAT=ppi; SUFFIX=; fi
OUT="pdisjoint_runs_v2/phosphorylation_st__${DS}__${CD}__split${SP}${SUFFIX}.json"
[ -f "$OUT" ] && { echo "already done: $OUT"; exit 0; }

echo "task $SLURM_ARRAY_TASK_ID -> phosphorylation_st / $DS / $CD / split $SP / $FEAT"
echo "Host: $(hostname)  Time: $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python -c "import torch; m = torch.nn.GRU(4, 4).cuda(); m(torch.zeros(1, 1, 4).cuda()); torch.cuda.synchronize()" \
    || { echo "GPU not supported by this torch build"; exit 1; }
python scripts/train_mouse.py --ptm phosphorylation_st --dataset "$DS" --cond "$CD" \
       --split-seed "$SP" --feat "$FEAT" --n-models 2 || { echo "training failed"; exit 1; }
[ -f "$OUT" ] || { echo "no result written"; exit 1; }
echo "Done at $(date)"
