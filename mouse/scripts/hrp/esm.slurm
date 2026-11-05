#!/bin/bash
### Mouse replication (PREREG_mouse.md): frozen mean-pooled ESM-2 650M vectors, repository extract_esm.py unchanged.
#SBATCH --job-name=mouse_esm
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --exclude=mantis-048,mantis-075
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jli89@tulane.edu
#SBATCH -e /home/FCAM/juli/HRP/mouse_v1/logs/esm_%j.err
#SBATCH -o /home/FCAM/juli/HRP/mouse_v1/logs/esm_%j.log
set -eo pipefail
export PTM_AUDIT_BASE=/home/FCAM/juli/HRP/mouse_v1
source /home/FCAM/juli/HRP/my-venv/bin/activate
cd "$PTM_AUDIT_BASE"
echo "Host: $(hostname)  Time: $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python -c "import torch; m = torch.nn.GRU(4, 4).cuda(); m(torch.zeros(1, 1, 4).cuda()); torch.cuda.synchronize()"
# extract_esm.py writes to $PTM_AUDIT_BASE/notebooks/, never to the human features
python /home/FCAM/juli/HRP/ptm-audit/extract_esm.py --batch-tokens 8192
echo "Done at $(date)"
