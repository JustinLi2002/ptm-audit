#!/usr/bin/env python3
"""
map_string_to_uniprot.py — the missing link in the embedding pipeline.

    node2vec_train.py        ->  ppi_node2vec.emb*.txt        (STRING ENSP ids)
    map_string_to_uniprot.py ->  node2vec_with_uniprot.csv    (UniProt accessions)  <- here
    build_protein_features.py->  protein_features_ppi.npy + protein_ids_ppi.json

This step was originally performed interactively and is reconstructed here. Because
several embedding files and several id-mapping tables exist in the working
directory, and because the handling of many-to-one mappings was not recorded, the
script has a --verify mode: it rebuilds the matrix under each duplicate-handling
policy and reports which one reproduces an existing protein_features_ppi.npy.
Run --verify first, then use the policy it identifies.

Usage:
    python map_string_to_uniprot.py --emb results/ppi_node2vec.emb.txt \\
        --idmap results/idmapping_2025_07_09.tsv \\
        --verify notebooks/protein_features_ppi.npy notebooks/protein_ids_ppi.json

    python map_string_to_uniprot.py --emb ... --idmap ... \\
        --policy first --out node2vec_with_uniprot.csv
"""

import argparse
import json
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd


def read_w2v(path):
    """word2vec text format: header 'n dim', then 'id v1 ... vdim' per line."""
    vecs, ids = [], []
    with open(path) as fh:
        head = fh.readline().split()
        if len(head) == 2 and head[0].isdigit():
            n, dim = int(head[0]), int(head[1])
        else:                                   # no header; rewind
            fh.seek(0)
            n, dim = None, None
        for line in fh:
            p = line.rstrip("\n").split()
            if len(p) < 3:
                continue
            ids.append(p[0])
            vecs.append(np.asarray(p[1:], dtype=np.float32))
    V = np.vstack(vecs)
    if dim and V.shape[1] != dim:
        print(f"  WARN: header says dim={dim}, parsed {V.shape[1]}", file=sys.stderr)
    print(f"  embedding: {V.shape[0]} nodes x {V.shape[1]} dims", file=sys.stderr)
    return ids, V


def strip_taxon(sid):
    """'9606.ENSP00000000233' -> 'ENSP00000000233'"""
    return sid.split(".", 1)[1] if "." in sid else sid


def read_idmap(path):
    """UniProt id-mapping export. Auto-detects the two relevant columns."""
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    cols = list(df.columns)
    # UniProt's export is usually 'From'/'Entry' or 'From'/'To'
    src = next((c for c in cols if c.lower() in ("from", "yourlist", "query")), cols[0])
    dst = next((c for c in cols if c.lower() in ("entry", "to", "accession")),
               cols[1] if len(cols) > 1 else cols[0])
    print(f"  idmap columns: '{src}' -> '{dst}'  ({len(df)} rows)", file=sys.stderr)
    m = defaultdict(list)
    for a, b in zip(df[src].values, df[dst].values):
        if a and b:
            m[strip_taxon(a.strip())].append(b.strip())
    return m


def build(ids, V, idmap, policy):
    """Return (accessions, matrix) under the given duplicate policy."""
    per = defaultdict(list)
    for sid, row in zip(ids, V):
        for acc in idmap.get(strip_taxon(sid), ()):
            per[acc].append(row)
    accs = sorted(per)
    if policy == "mean":
        M = np.vstack([np.mean(per[a], axis=0) for a in accs])
    elif policy == "first":
        M = np.vstack([per[a][0] for a in accs])
    elif policy == "unique":                     # drop accessions with >1 source node
        accs = [a for a in accs if len(per[a]) == 1]
        M = np.vstack([per[a][0] for a in accs])
    else:
        raise ValueError(policy)
    return accs, M.astype(np.float32)


def compare(accs, M, ref_npy, ref_ids):
    """Report how closely a rebuild matches an existing feature matrix."""
    R = np.load(ref_npy)
    with open(ref_ids) as fh:
        rids = json.load(fh)
    same_n = (len(accs) == len(rids))
    same_order = same_n and accs == list(rids)
    inter = set(accs) & set(rids)
    out = [f"n={len(accs)} vs ref {len(rids)}", f"shared accessions {len(inter)}"]
    if inter:
        ia = {a: i for i, a in enumerate(accs)}
        ir = {a: i for i, a in enumerate(rids)}
        sub = sorted(inter)[:2000]
        d = np.abs(M[[ia[a] for a in sub]] - R[[ir[a] for a in sub]])
        out.append(f"max|diff| on {len(sub)} shared rows = {d.max():.2e}")
        out.append(f"identical rows: {int((d.max(axis=1) < 1e-6).sum())}/{len(sub)}")
    out.append(f"same length: {same_n}; same order: {same_order}")
    return "  " + "\n  ".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb", required=True, help="word2vec-format embedding")
    ap.add_argument("--idmap", required=True, help="UniProt id-mapping export (.tsv)")
    ap.add_argument("--policy", default="mean", choices=["mean", "first", "unique"],
                    help="how to handle several STRING nodes mapping to one accession")
    ap.add_argument("--out", help="output CSV (omit to only report)")
    ap.add_argument("--verify", nargs=2, metavar=("NPY", "IDS_JSON"),
                    help="compare rebuilds against an existing feature matrix")
    a = ap.parse_args()

    for f in (a.emb, a.idmap):
        if not os.path.exists(f):
            sys.exit(f"missing: {f}")

    print(f"reading {a.emb}", file=sys.stderr)
    ids, V = read_w2v(a.emb)
    print(f"reading {a.idmap}", file=sys.stderr)
    idmap = read_idmap(a.idmap)

    mapped = sum(1 for s in ids if strip_taxon(s) in idmap)
    print(f"  {mapped}/{len(ids)} embedding nodes have a UniProt mapping "
          f"({100*mapped/len(ids):.1f}%)\n", file=sys.stderr)

    if a.verify:
        for pol in ("mean", "first", "unique"):
            accs, M = build(ids, V, idmap, pol)
            print(f"policy = {pol}")
            print(compare(accs, M, a.verify[0], a.verify[1]))
            print()
        print("Use the policy whose rebuild matches in length, order and values.\n"
              "If none matches exactly, the archived matrix should be treated as the\n"
              "canonical artefact and deposited alongside the code.")
        return

    accs, M = build(ids, V, idmap, a.policy)
    print(f"built {M.shape[0]} accessions x {M.shape[1]} dims "
          f"(policy={a.policy})", file=sys.stderr)
    if a.out:
        df = pd.DataFrame(M, index=accs)
        df.index.name = "uniprot"
        df.to_csv(a.out)
        print(f"wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
