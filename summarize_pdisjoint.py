import json, glob, collections
import numpy as np
R = collections.defaultdict(dict)
for f in glob.glob('pdisjoint_runs/*.json'):
    d = json.load(open(f))
    R[(d['ptm'], d['dataset'], d['cond'])][d['split_seed']] = d
print(f"完成 {len(glob.glob('pdisjoint_runs/*.json'))}/96\n")
ptms = sorted({k[0] for k in R})
print(f"{'PTM':20s} {'条件':20s} {'AUROC (跨split)':>22s} {'AUPRC':>8s} {'n':>3s}")
print('-'*78)
for p in ptms:
    for ds in ['replica','rebuilt']:
        for c in ['baseline','ppi']:
            runs = R.get((p,ds,c))
            if not runs: continue
            a = [r['auroc'] for r in runs.values()]
            u = [r['auprc'] for r in runs.values()]
            print(f"{p:20s} {ds+'/'+c:20s} {np.mean(a):9.4f} ± {np.std(a):.4f}"
                  f"  {np.mean(u):8.4f} {len(a):3d}")
    # PPI 增益
    for ds in ['replica','rebuilt']:
        b, q = R.get((p,ds,'baseline')), R.get((p,ds,'ppi'))
        if b and q:
            common = set(b) & set(q)
            if common:
                g = [q[s]['auroc']-b[s]['auroc'] for s in common]
                print(f"{'':20s} {ds+' PPI增益':20s} {np.mean(g):+9.4f} ± {np.std(g):.4f}"
                      f"{'':10s} {len(g):3d}")
    print()
