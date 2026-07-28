#!/bin/bash
### baseline / +PPI / permuted-PPI on the authors' final released dataset.
### 8 PTMs x 3 conditions x 2 validation seeds = 48 cells, 2 models each.
#SBATCH --job-name=alldata
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=98G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jli89@tulane.edu
#SBATCH --time=23:00:00
#SBATCH --array=0-47%8
#SBATCH -e logs/alldata_%A_%a.err
#SBATCH -o logs/alldata_%A_%a.log

source /home/FCAM/juli/HRP/my-venv/bin/activate
cd /home/FCAM/juli/HRP
mkdir -p logs alldata_runs

PTMS=(acetylation_k glycosylation_n methylation_k methylation_r \
      phosphorylation_st phosphorylation_y sumoylation_k ubiquitination_k)
CONDS=(baseline ppi shuffled)
VSEEDS=(0 1)

i=$SLURM_ARRAY_TASK_ID
v=$(( i % 2 )); i=$(( i / 2 ))
c=$(( i % 3 )); i=$(( i / 3 ))
p=$(( i % 8 ))

PTM=${PTMS[$p]}; CD=${CONDS[$c]}; VS=${VSEEDS[$v]}
OUT="alldata_runs/${PTM}__${CD}__val${VS}__f90.json"
if [ -f "$OUT" ]; then echo "already done: $OUT"; exit 0; fi

echo "task $SLURM_ARRAY_TASK_ID -> $PTM / $CD / val seed $VS"
echo "Host: $(hostname)  Time: $(date)"

python train_alldata.py --ptm "$PTM" --cond "$CD" --val-seed "$VS" --n-models 2

echo "Done at $(date)"
