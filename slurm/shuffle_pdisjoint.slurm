#!/bin/bash
### Shuffle control under the protein-disjoint split.
### 8 PTMs x 3 split seeds = 24 cells, rebuilt training data only.
### Asks whether the +0.021 AUROC / +0.050 AUPRC gain seen with corrected
### sampling survives permutation of the protein-to-embedding mapping.
#SBATCH --job-name=pdisj_shuf
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=98G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=justin.2.li@uconn.edu
#SBATCH --time=12:00:00
#SBATCH --array=0-23%8
#SBATCH -e logs/pdisj_shuf_%A_%a.err
#SBATCH -o logs/pdisj_shuf_%A_%a.log

source /home/FCAM/juli/HRP/my-venv/bin/activate
cd /home/FCAM/juli/HRP
mkdir -p logs pdisjoint_runs

PTMS=(acetylation_k glycosylation_n methylation_k methylation_r \
      phosphorylation_st phosphorylation_y sumoylation_k ubiquitination_k)
SPLITS=(0 1 2)

i=$SLURM_ARRAY_TASK_ID
s=$(( i % 3 )); i=$(( i / 3 ))
p=$(( i % 8 ))

PTM=${PTMS[$p]}; SP=${SPLITS[$s]}
OUT="pdisjoint_runs/${PTM}__rebuilt__shuffled__split${SP}.json"

if [ -f "$OUT" ]; then
    echo "already done: $OUT"; exit 0
fi

echo "task $SLURM_ARRAY_TASK_ID -> $PTM / rebuilt / shuffled / split $SP"
echo "Host: $(hostname)  Time: $(date)"

python train_pdisjoint.py --ptm "$PTM" --dataset rebuilt --cond shuffled \
       --split-seed "$SP" --n-models 2

echo "Done at $(date)"
