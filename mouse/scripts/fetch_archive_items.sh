#!/bin/bash
# One connection (one Duo push): the mouse_v1 files that exist only on HRP and are needed for the v1.6 archive.
# - matched/ (ratio-matched training and validation sets per split, matched_summary.tsv, manifest.sha256)
# - notebooks/ (ESM-2 feature matrix and protein ids)
# - checkpoints_pdisjoint/ (model weights written by train_pdisjoint.py), with sha256 computed on HRP
# plus an inventory of mouse_v1. Run in the foreground from mouse_replication_20261102/.
set -u
mkdir -p fetch_archive
ssh -n USER@HPC_LOGIN_HOST 'cd /home/FCAM/juli/HRP/mouse_v1 &&
  { echo "== du"; du -sh * 2>/dev/null; echo "== checkpoints"; ls checkpoints_pdisjoint 2>/dev/null | wc -l;
    echo "== matched"; find matched -type f | sort; } > archive_inventory.txt;
  ( cd checkpoints_pdisjoint 2>/dev/null && sha256sum *.pt ) > checkpoints.sha256;
  ( cd matched && find . -name "*.tsv" -path "./split*" | sort | xargs sha256sum ) > matched_files.sha256;
  tar czf - archive_inventory.txt checkpoints.sha256 matched_files.sha256 matched notebooks \
    checkpoints_pdisjoint' > fetch_archive/mouse_archive_items.tgz
tar --force-local -xzf fetch_archive/mouse_archive_items.tgz -C fetch_archive
cat fetch_archive/archive_inventory.txt
echo "== local check"
( cd fetch_archive/checkpoints_pdisjoint && sha256sum -c --quiet ../checkpoints.sha256 && echo "checkpoints: sha256 OK" )
( cd fetch_archive/matched && sha256sum -c --quiet ../matched_files.sha256 && echo "matched: sha256 OK" )
( cd fetch_archive && sha256sum -c --quiet ../fetch/notebooks.sha256 && echo "notebooks: sha256 OK (matches the hash taken at the results fetch)" )
