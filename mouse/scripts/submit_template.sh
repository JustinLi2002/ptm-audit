#!/bin/bash
# Mouse replication (PREREG_mouse.md): one login-node connection that unpacks mouse_v1 and submits
# prep -> (esm) -> 27-cell training array. Run as: ssh USER@HPC_LOGIN_HOST 'bash -s' < submit_mouse_v1.sh
set -eo pipefail
cd /home/FCAM/juli/HRP
if [ -e mouse_v1 ]; then
    # a previous connection got this far: report its state instead of submitting again
    echo "mouse_v1 already exists; nothing submitted"
    cat mouse_v1/jobs.txt 2>/dev/null || echo "no jobs.txt"
    squeue -u juli -o "%.12i %.14j %.8T %.10M %.20R" | head -40
    ls -la mouse_v1 mouse_v1/logs 2>/dev/null
    exit 0
fi
for f in ptm-audit/train_pdisjoint.py ptm-audit/make_pdisjoint_split.py ptm-audit/extract_esm.py \
         /home/FCAM/juli/anaconda3/envs/deepmvp/bin/cd-hit my-venv/bin/activate; do
    [ -e "$f" ] || { echo "missing: $f"; exit 1; }
done
base64 -d > mouse_v1_bundle.tgz <<'B64'
@@B64@@
B64
echo "@@SHA@@  mouse_v1_bundle.tgz" | sha256sum -c -
tar xzf mouse_v1_bundle.tgz
mkdir -p mouse_v1/logs
echo "== repository scripts used (sha256, LF)"
sha256sum ptm-audit/train_pdisjoint.py ptm-audit/make_pdisjoint_split.py ptm-audit/extract_esm.py
source my-venv/bin/activate
PTM_AUDIT_BASE=/home/FCAM/juli/HRP/mouse_v1 python mouse_v1/scripts/train_mouse.py --check
ls ~/.cache/torch/hub/checkpoints/ 2>/dev/null | grep -i "esm2_t33_650M" || echo "WARN: ESM-2 650M weights not in the torch hub cache"
P=$(sbatch --parsable mouse_v1/scripts/prep.slurm)
E=$(sbatch --parsable mouse_v1/scripts/esm.slurm)
T=$(sbatch --parsable --dependency=afterok:$P:$E mouse_v1/scripts/train.slurm)
echo "jobs: prep $P  esm $E  train $T"
echo "$P $E $T" > mouse_v1/jobs.txt
squeue -u juli -o "%.12i %.14j %.8T %.10M %.20R" | head -20
