#!/bin/bash
# One connection (one Duo push): job states, matched-run outputs, logs, and the existing T/U run summaries of the two
# tasks, into fetch/. Run in the foreground from this directory, after the END mail for array 4406684.
set -eu
mkdir -p fetch
ssh -n USER@HPC_LOGIN_HOST 'cd /home/FCAM/juli/HRP &&
  sacct -j 4406684,4406841 -n -P --format=JobID,NodeList,State,Elapsed,ExitCode > ratio_matched_v4/sacct.txt;
  tar czf - ratio_matched_v4/sacct.txt matched_v4/matched_summary.tsv matched_v4/manifest.sha256 \
    pdisjoint_runs_matched_v4 logs/matched4_4406684_* logs/matched4r_4406841_* \
    pdisjoint_runs_v2/acetylation_k__*.json pdisjoint_runs_v2/phosphorylation_st__*.json' > fetch/fetch_real.tgz
tar --force-local -xzf fetch/fetch_real.tgz -C fetch
cat fetch/ratio_matched_v4/sacct.txt
