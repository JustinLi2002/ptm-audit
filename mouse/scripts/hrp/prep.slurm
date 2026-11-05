#!/bin/bash
### Mouse replication (PREREG_mouse.md): CD-HIT clusters, three protein-disjoint splits, matched construction.
#SBATCH --job-name=mouse_prep
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jli89@tulane.edu
#SBATCH -e /home/FCAM/juli/HRP/mouse_v1/logs/prep_%j.err
#SBATCH -o /home/FCAM/juli/HRP/mouse_v1/logs/prep_%j.log
set -eo pipefail
export PTM_AUDIT_BASE=/home/FCAM/juli/HRP/mouse_v1
cd "$PTM_AUDIT_BASE"
echo "Host: $(hostname)  Time: $(date)"

# same CD-HIT call as slurm/cdhit.slurm (-c 0.4 -n 2), over every protein of the unrestricted construction;
# the conda env's libraries stand in for `conda activate deepmvp`
export LD_LIBRARY_PATH=/home/FCAM/juli/anaconda3/envs/deepmvp/lib:${LD_LIBRARY_PATH:-}
/home/FCAM/juli/anaconda3/envs/deepmvp/bin/cd-hit -i pdisjoint/proteins.fasta -o pdisjoint/clust40 \
    -c 0.4 -n 2 -g 1 -d 0 -M 40000 -T 16

source /home/FCAM/juli/HRP/my-venv/bin/activate
for s in 0 1 2; do
    python /home/FCAM/juli/HRP/ptm-audit/make_pdisjoint_split.py --clstr pdisjoint/clust40.clstr \
        --sites pdisjoint/sites_for_split.csv --out pdisjoint/split_seed$s.csv --seed $s
done
python scripts/build_matched_mouse.py
sha256sum pdisjoint/split_seed*.csv pdisjoint/clust40.clstr > pdisjoint/manifest.sha256
echo "Done at $(date)"
