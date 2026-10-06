#!/usr/bin/env python3
"""Rerun only FN sensitivity cells affected by deduplicating methylation.

Uses the archived Monte Carlo algorithm, grid, 200 draws, six bins, and seed 0.
Unchanged tasks are copied from the published output. alpha=1 must reproduce
the corresponding original cell before any revised result can be assembled.
"""
import argparse
import importlib.util
import sys
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--shard',type=int,required=True)
    a=ap.parse_args()
    arm=('ppi','esm')[a.shard//6]
    task=('methylation_k','methylation_r')[(a.shard%6)//3]
    split=a.shard%3
    archive=pd.read_csv(a.root/f'ptm-audit/results/fn_sensitivity_{arm}.tsv',sep='\t')
    alphas=sorted(archive.loc[archive['mode'].eq('alpha'),'alpha'].unique())
    targets=sorted(archive.loc[archive['mode'].eq('inflation'),'target_inflation'].unique())
    script=a.root/'ptm-audit/analysis/fn_sensitivity.py'
    spec=importlib.util.spec_from_file_location('fn_original',script)
    fn=importlib.util.module_from_spec(spec);spec.loader.exec_module(fn)
    fn.load_project_modules(arm)
    # The first file already contains the full K/R set. Verify and take a union.
    frames=[pd.read_csv(a.root/f'rebuilt/{t}_all.tsv',sep='\t',usecols=['protein','pos','y'])
            for t in ('methylation_k','methylation_r')]
    positives=[d[d.y.eq(1)][['protein','pos']].drop_duplicates() for d in frames]
    merged=pd.concat(positives).drop_duplicates(['protein','pos']).groupby('protein').size().to_dict()
    fn.depths=lambda: ({},{'methyl':merged})
    fn.PTMS=[task];fn.SPLITS=(split,)
    dest=a.out/f'shard{a.shard:02d}'
    df=fn.run(arm,alphas,200,6,0,dest,targets)
    old=archive[(archive.task==task)&(archive['split']==split)&(archive['mode']=='alpha')&(archive.alpha==1)].iloc[0]
    new=df[(df['mode']=='alpha')&(df.alpha==1)].iloc[0]
    for metric in ('d_auroc','d_auprc'):
        assert np.isclose(old[metric],new[metric],atol=1e-12,rtol=0),(task,split,metric,old[metric],new[metric])
    (dest/'VALIDATED').write_text('alpha=1 reproduces archived AUROC and AUPRC within 1e-12\n')
    print('SENSITIVITY_SHARD_VALIDATED',a.shard,arm,task,split,flush=True)


if __name__=='__main__': main()
