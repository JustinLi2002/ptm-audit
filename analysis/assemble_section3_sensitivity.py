from pathlib import Path
import pandas as pd
root=Path('/home/FCAM/juli/HRP')
work=root/'section3_revision_20260920'
out=work/'results'
for arm in ('ppi','esm'):
    old=pd.read_csv(root/'ptm-audit/results'/f'fn_sensitivity_{arm}.tsv',sep='\t')
    pieces=[]
    for shard in range(12):
        d=work/'fn_shards'/f'shard{shard:02}'
        assert (d/'VALIDATED').exists(),d
        f=d/f'fn_sensitivity_{arm}.tsv'
        if f.exists(): pieces.append(pd.read_csv(f,sep='\t'))
    new=pd.concat(pieces,ignore_index=True)
    assert len(new)==len(old[old.task.isin(['methylation_k','methylation_r'])])
    keep=old[~old.task.isin(['methylation_k','methylation_r'])]
    combined=pd.concat([keep,new],ignore_index=True).sort_values(['task','split','mode','alpha','target_inflation'])
    combined.to_csv(out/f'fn_sensitivity_{arm}.tsv',sep='\t',index=False)
    print(arm,len(old),len(combined))
