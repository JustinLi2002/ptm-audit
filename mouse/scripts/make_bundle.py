#!/usr/bin/env python3
"""Assemble the mouse_v1 tree for HRP from build/ and scripts/hrp/, tar it, and write the one-connection submit
script (base64 tarball embedded; checked by sha256 on arrival).

Usage: python make_bundle.py  (run from mouse_replication_20261102/)"""
import base64
import hashlib
import io
import shutil
import tarfile
from pathlib import Path

import pandas as pd

ROOT = Path.cwd()
B = ROOT / "build"
TREE = B / "hrp" / "mouse_v1"
if TREE.exists():
    shutil.rmtree(TREE)
for sub in ("replica", "rebuilt", "pdisjoint", "scripts", "notes"):
    (TREE / sub).mkdir(parents=True)

for c in ("replica", "rebuilt"):
    # pandas on Windows writes CRLF; the HRP copies are LF
    (TREE / c / "phosphorylation_st_all.tsv").write_bytes(
        (B / c / "phosphorylation_st_all.tsv").read_bytes().replace(b"\r\n", b"\n"))
U = pd.read_csv(B / "rebuilt" / "phosphorylation_st_all.tsv", sep="\t")
T = pd.read_csv(B / "replica" / "phosphorylation_st_all.tsv", sep="\t")
pos = lambda d: set(zip(d.protein[d.y == 1], d.pos[d.y == 1]))
assert pos(T) == pos(U)                                   # the constructions share their positives
assert set(T.protein) == set(U.protein) == set(U[U.y == 1].protein)
prots = sorted(U.protein.unique())

seqs, acc = {}, None
for line in open(B / "reference.fasta", encoding="utf-8"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        acc = line.split("|")[1]
        seqs[acc] = []
    else:
        seqs[acc].append(line)
with open(TREE / "pdisjoint" / "proteins.fasta", "w", encoding="utf-8", newline="\n") as fh:
    for p in prots:
        fh.write(f">{p}\n" + "\n".join(seqs[p]) + "\n")
U.rename(columns={"protein": "accession", "y": "label"}).assign(ptm_type="phosphorylation_st")[
    ["accession", "ptm_type", "label"]].to_csv(TREE / "pdisjoint" / "sites_for_split.csv", index=False,
                                               lineterminator="\n")

for f in (ROOT / "scripts" / "hrp").iterdir():
    shutil.copyfile(f, TREE / "scripts" / f.name)
for f in (ROOT / "PREREG_mouse.md", B / "construction_prediction.txt", B / "mapping_report.json"):
    (TREE / "notes" / f.name).write_bytes(f.read_bytes().replace(b"\r\n", b"\n"))
for f in TREE.rglob("*"):
    if f.is_file():
        assert b"\r" not in f.read_bytes(), f

buf = io.BytesIO()
with tarfile.open(fileobj=buf, mode="w:gz") as tar:
    tar.add(TREE, arcname="mouse_v1")
blob = buf.getvalue()
sha = hashlib.sha256(blob).hexdigest()
(B / "hrp" / "mouse_v1_bundle.tgz").write_bytes(blob)
b64 = base64.encodebytes(blob).decode()
script = (ROOT / "scripts" / "submit_template.sh").read_text(encoding="utf-8")
script = script.replace("@@SHA@@", sha).replace("@@B64@@\n", b64)
(B / "hrp" / "submit_mouse_v1.sh").write_text(script, encoding="utf-8", newline="\n")
print(f"{len(prots)} proteins; tarball {len(blob) / 1e6:.2f} MB, base64 {len(b64) / 1e6:.2f} MB, sha256 {sha}")
