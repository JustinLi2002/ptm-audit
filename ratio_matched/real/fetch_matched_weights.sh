#!/bin/bash
# One connection (one Duo push): the human ratio-matched model weights for the v1.6 data record, with sha256 computed
# on HRP, plus the repository commit they were trained from. Run in the foreground from archive_v1.6/.
set -u
mkdir -p fetch_hrp
ssh -n USER@HPC_LOGIN_HOST 'cd /home/FCAM/juli/HRP &&
  { echo "== repo HEAD"; git -C ptm-audit rev-parse HEAD; echo "== checkpoints"; ls checkpoints_matched_v4 | wc -l;
    du -sh checkpoints_matched_v4; echo "== train script"; sha256sum ptm-audit/train_pdisjoint.py ratio_matched_v4/train_matched.py; } > /home/FCAM/juli/HRP/matched_v4/archive_inventory.txt 2>&1;
  ( cd checkpoints_matched_v4 && sha256sum *.pt ) > /home/FCAM/juli/HRP/matched_v4/checkpoints_matched.sha256;
  tar czf - -C /home/FCAM/juli/HRP/matched_v4 archive_inventory.txt checkpoints_matched.sha256 -C /home/FCAM/juli/HRP checkpoints_matched_v4' > fetch_hrp/matched_weights.tgz
ls -la fetch_hrp/matched_weights.tgz
tar --force-local -xzf fetch_hrp/matched_weights.tgz -C fetch_hrp
cat fetch_hrp/archive_inventory.txt
echo "== local check"
( cd fetch_hrp/checkpoints_matched_v4 && sha256sum -c --quiet ../checkpoints_matched.sha256 && echo "checkpoints: sha256 OK" && ls | wc -l )
