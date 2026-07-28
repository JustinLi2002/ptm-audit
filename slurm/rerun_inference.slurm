#!/bin/bash
#SBATCH --job-name=deepmvp_ppi
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=60G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jli89@tulane.edu
#SBATCH -e inference.err
#SBATCH -o inference.log

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

source /home/FCAM/juli/HRP/my-venv/bin/activate

echo "=== GPU INFO ==="
nvidia-smi
echo "================"

cd /home/FCAM/juli/HRP
python rerun_inference.py --cond all

echo "Done at" `date`
