#!/usr/bin/env python3
"""Mean-pooled ESM-2 650M embeddings, one vector per protein.

A frozen mean-pooled protein language model embedding is constant across all
candidate sites within a protein -- structurally identical to the interaction
embedding analysed in the main experiments. This tests whether the mechanism
is a property of hand-crafted protein-level features or of protein-constant
inputs in general.

Sequences longer than the 1022-residue context are embedded in overlapping
chunks and pooled with weights proportional to the number of residues each
chunk contributes, which reduces to a plain mean for short sequences.

Output matches the interaction feature format so train_pdisjoint.py can use
either: protein_features_esm.npy (n x 1280) + protein_ids_esm.json.
"""
import argparse
import json
import os

import numpy as np
import torch

BASE = os.environ.get("PTM_AUDIT_BASE", "/home/FCAM/juli/HRP")
FASTA = f"{BASE}/pdisjoint/proteins.fasta"
OUT_NPY = f"{BASE}/notebooks/protein_features_esm.npy"
OUT_IDS = f"{BASE}/notebooks/protein_ids_esm.json"
WINDOW, OVERLAP = 1022, 256


def read_fasta(path):
    ids, seqs, cur, buf = [], [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            if cur is not None:
                ids.append(cur); seqs.append(''.join(buf))
            cur, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if cur is not None:
        ids.append(cur); seqs.append(''.join(buf))
    return ids, seqs


def chunks(seq):
    """Overlapping windows covering the whole sequence."""
    if len(seq) <= WINDOW:
        return [(seq, len(seq))]
    step = WINDOW - OVERLAP
    out, start = [], 0
    while start < len(seq):
        piece = seq[start:start + WINDOW]
        out.append((piece, len(piece)))
        if start + WINDOW >= len(seq):
            break
        start += step
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--batch-tokens', type=int, default=8192)
    ap.add_argument('--limit', type=int, default=0, help='debug: first N proteins')
    a = ap.parse_args()

    dev = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'device={dev}', flush=True)
    model, alphabet = torch.hub.load_state_dict_from_url, None
    import esm
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(dev).eval()
    bc = alphabet.get_batch_converter()
    LAYER, DIM = 33, 1280

    ids, seqs = read_fasta(FASTA)
    if a.limit:
        ids, seqs = ids[:a.limit], seqs[:a.limit]
    print(f'{len(ids)} sequences, max len {max(len(s) for s in seqs)}', flush=True)

    # flatten to chunks, sort by length so batches are homogeneous
    flat = [(i, piece, n) for i, s in enumerate(seqs) for piece, n in chunks(s)]
    flat.sort(key=lambda t: len(t[1]))
    print(f'{len(flat)} chunks', flush=True)

    acc = np.zeros((len(ids), DIM), dtype=np.float64)
    wsum = np.zeros(len(ids), dtype=np.float64)

    batch, ntok, done = [], 0, 0
    def flush(batch):
        if not batch:
            return
        _, _, toks = bc([(f'p{k}', p) for k, (i, p, n) in enumerate(batch)])
        with torch.no_grad():
            rep = model(toks.to(dev), repr_layers=[LAYER])['representations'][LAYER]
        for k, (i, p, n) in enumerate(batch):
            v = rep[k, 1:len(p) + 1].mean(0).float().cpu().numpy()
            acc[i] += v * n
            wsum[i] += n

    for item in flat:
        L = len(item[1]) + 2
        if batch and ntok + L > a.batch_tokens:
            flush(batch); done += len(batch)
            if done % 2000 < len(batch):
                print(f'  {done}/{len(flat)} chunks', flush=True)
            batch, ntok = [], 0
        batch.append(item); ntok += L
    flush(batch); done += len(batch)
    print(f'  {done}/{len(flat)} chunks', flush=True)

    assert (wsum > 0).all(), 'some protein got no chunk'
    feats = (acc / wsum[:, None]).astype(np.float32)
    os.makedirs(os.path.dirname(OUT_NPY), exist_ok=True)
    np.save(OUT_NPY, feats)
    json.dump(ids, open(OUT_IDS, 'w'))
    print(f'wrote {feats.shape} to {OUT_NPY}', flush=True)

    # sanity: no collapse
    rng = np.random.default_rng(0)
    a_, b_ = rng.integers(0, len(ids), 5000), rng.integers(0, len(ids), 5000)
    cs = (feats[a_] * feats[b_]).sum(1) / (
        np.linalg.norm(feats[a_], axis=1) * np.linalg.norm(feats[b_], axis=1))
    print(f'mean cosine of random pairs = {cs.mean():.3f} (collapse would be ~1)')


if __name__ == '__main__':
    main()
