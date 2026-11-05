#!/usr/bin/env python3
"""Mouse replication (PREREG_mouse.md): runs the repository's train_pdisjoint.py on the mouse_v1 tree.

The repository script is not changed. Two strings are replaced, so that 'matched' is accepted and its training
file is read from matched/split{s}/.

All paths come from PTM_AUDIT_BASE, which must be the mouse_v1 directory. Otherwise the script would read and
overwrite the human runs.

    python train_mouse.py --check     # verify the guard and the patch, then exit
"""
import os
import sys
from pathlib import Path

BASE = os.environ.get("PTM_AUDIT_BASE", "")
if Path(BASE).name != "mouse_v1":
    sys.exit(f"PTM_AUDIT_BASE must point to mouse_v1, got {BASE!r}")
SRC = Path(__file__).resolve().parents[3] / "train_pdisjoint.py"
src = SRC.read_text()
EDITS = [
    ('choices=["replica", "rebuilt"]', 'choices=["replica", "rebuilt", "matched"]'),
    ('train_src = pd.read_csv(f"{BASE}/{a.dataset}/{a.ptm}_all.tsv", sep="\\t")',
     'train_src = pd.read_csv(f"{BASE}/matched/split{a.split_seed}/{a.ptm}_all.tsv" '
     'if a.dataset == "matched" else f"{BASE}/{a.dataset}/{a.ptm}_all.tsv", sep="\\t")'),
]
for old, new in EDITS:
    assert src.count(old) == 1, old
    src = src.replace(old, new)
assert 'BASE = os.environ.get("PTM_AUDIT_BASE"' in src

if "--check" in sys.argv:
    print("patch applies:", len(EDITS), "edits; base", BASE)
    sys.exit(0)
sys.argv[0] = str(SRC)
exec(compile(src, str(SRC), "exec"), {"__name__": "__main__", "__file__": str(SRC)})
