#!/bin/bash
#SBATCH --job-name=esmfeat
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jli89@tulane.edu
#SBATCH -e logs/esmfeat_%j.err
#SBATCH -o logs/esmfeat_%j.log

source /home/FCAM/juli/HRP/my-venv/bin/activate
cd /home/FCAM/juli/HRP
mkdir -p logs
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python ptm-audit/extract_esm.py --batch-tokens 8192

