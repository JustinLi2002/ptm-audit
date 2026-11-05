#!/usr/bin/env python3
"""train_pdisjoint.py with a third training construction, 'matched' (PREREG_real.md).

Runs the repository script with four strings replaced: the training-data path for 'matched', the allowed
--dataset values and the two output directories. The repository script itself is not changed.

    python train_matched.py --check            # verify the patch applies, then exit
"""
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[2] / "train_pdisjoint.py"
src = SRC.read_text()
EDITS = [
    ('OUTDIR = f"{BASE}/pdisjoint_runs_v2"', 'OUTDIR = f"{BASE}/pdisjoint_runs_matched_v4"'),
    ('CKPT_DIR = f"{BASE}/checkpoints_pdisjoint"', 'CKPT_DIR = f"{BASE}/checkpoints_matched_v4"'),
    ('choices=["replica", "rebuilt"]', 'choices=["replica", "rebuilt", "matched"]'),
    ('train_src = pd.read_csv(f"{BASE}/{a.dataset}/{a.ptm}_all.tsv", sep="\\t")',
     'train_src = pd.read_csv(f"{BASE}/matched_v4/split{a.split_seed}/{a.ptm}_all.tsv" '
     'if a.dataset == "matched" else f"{BASE}/{a.dataset}/{a.ptm}_all.tsv", sep="\\t")'),
]
for old, new in EDITS:
    assert src.count(old) == 1, old
    src = src.replace(old, new)

if "--check" in sys.argv:
    print("patch applies:", len(EDITS), "edits")
    sys.exit(0)
sys.argv[0] = str(SRC)
exec(compile(src, str(SRC), "exec"), {"__name__": "__main__", "__file__": str(SRC)})
