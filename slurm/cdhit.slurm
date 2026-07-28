#!/bin/bash
#SBATCH --job-name=cdhit40
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jli89@tulane.edu
#SBATCH -e cdhit.err
#SBATCH -o cdhit.log

echo "Host: $(hostname)  Time: $(date)"
source ~/.bashrc
conda activate deepmvp
cd /home/FCAM/juli/HRP

/home/FCAM/juli/anaconda3/envs/deepmvp/bin/cd-hit -i pdisjoint/proteins.fasta -o pdisjoint/clust40 \
    -c 0.4 -n 2 -g 1 -d 0 -M 40000 -T 16

echo "Done at $(date)"


