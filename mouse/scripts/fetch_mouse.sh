#!/bin/bash
# One connection (one Duo push): job states, splits, matched summary, feature hashes, runs and logs of mouse_v1,
# into fetch/. Run in the foreground from mouse_replication_20261102/, after the END mail for the training array.
set -u
mkdir -p fetch
ssh -n USER@HPC_LOGIN_HOST 'cd /home/FCAM/juli/HRP/mouse_v1 &&
  squeue -u juli -o "%.14i %.14j %.8T %.10M %.20R" > squeue.txt;
  sacct -j $(cat jobs.txt | tr " " ",") -n -P --format=JobID,JobName,NodeList,State,Elapsed,ExitCode > sacct.txt;
  sha256sum notebooks/* > notebooks.sha256 2>/dev/null;
  tar czf - jobs.txt squeue.txt sacct.txt notebooks.sha256 pdisjoint/split_seed0.csv pdisjoint/split_seed1.csv \
    pdisjoint/split_seed2.csv pdisjoint/manifest.sha256 pdisjoint/clust40.clstr matched/matched_summary.tsv \
    matched/manifest.sha256 pdisjoint_runs_v2 logs' > fetch/fetch_mouse.tgz
tar --force-local -xzf fetch/fetch_mouse.tgz -C fetch
echo "== sacct"; cat fetch/sacct.txt
echo "== squeue"; cat fetch/squeue.txt
echo "== runs"; ls fetch/pdisjoint_runs_v2 2>/dev/null | grep -c "json$"
