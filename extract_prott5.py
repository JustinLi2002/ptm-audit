#!/usr/bin/env python3
"""Mean-pooled ProtT5-XL-U50 embeddings, one vector per protein.

A second protein language model, to test whether the reversal reported for
frozen ESM-2 embeddings is a property of protein-constant inputs in general
or of one model's representation.

The chunking and the length-weighted pooling are identical to extract_esm.py
(WINDOW 1022, OVERLAP 256, weights proportional to the residues each chunk
contributes).  ProtT5 uses relative position embeddings and has no hard
context limit, so the window is not required by the model; it is imposed so
that the two feature families differ only in the encoder.

ProtT5 input differs from ESM in three ways that fail silently if missed:
residues must be space-separated, the rare codes U/Z/O/B must be mapped to X,
and the tokenizer appends </s> but prepends nothing, so residue j sits at
position j rather than j+1.

Output matches the interaction and ESM feature format so train_pdisjoint.py
can use any of them: protein_features_prott5.npy (n x 1024) +
protein_ids_prott5.json.
"""
import argparse
import json
import os
import re

import numpy as np
import torch

BASE = os.environ.get("PTM_AUDIT_BASE", "/home/FCAM/juli/HRP")
FASTA = f"{BASE}/pdisjoint/proteins.fasta"
MODEL = "Rostlab/prot_t5_xl_half_uniref50-enc"
WINDOW, OVERLAP = 1022, 256
DIM = 1024


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
    """Overlapping windows covering the whole sequence. Same rule as ESM."""
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


def spaced(seq):
    """ProtT5 expects space-separated residues with rare codes mapped to X."""
    return ' '.join(re.sub(r'[UZOB]', 'X', seq))


def load_model(dev):
    from transformers import T5EncoderModel
    try:
        from transformers import T5Tokenizer as Tok
        tok = Tok.from_pretrained(MODEL, do_lower_case=False, legacy=True)
    except Exception:
        from transformers import AutoTokenizer as Tok
        tok = Tok.from_pretrained(MODEL, do_lower_case=False)
    try:
        model = T5EncoderModel.from_pretrained(MODEL, dtype=torch.float16)
    except TypeError:                       # transformers < 4.56
        model = T5EncoderModel.from_pretrained(MODEL, torch_dtype=torch.float16)
    if dev.type != 'cuda':                  # fp16 on CPU is not supported
        model = model.float()
    return tok, model.to(dev).eval()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--batch-tokens', type=int, default=4096,
                    help='ProtT5-XL is 3B parameters; activations are larger '
                         'than ESM-650M at the same token count')
    ap.add_argument('--limit', type=int, default=0,
                    help='debug: first N proteins, written to a _debug file')
    a = ap.parse_args()

    suffix = '_debug' if a.limit else ''
    out_npy = f"{BASE}/notebooks/protein_features_prott5{suffix}.npy"
    out_ids = f"{BASE}/notebooks/protein_ids_prott5{suffix}.json"

    dev = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'device={dev}', flush=True)
    if dev.type == 'cuda':
        print(f'gpu={torch.cuda.get_device_name(0)} '
              f'arch_list={torch.cuda.get_arch_list()}', flush=True)

    tok, model = load_model(dev)
    print(f'loaded {MODEL}', flush=True)

    ids, seqs = read_fasta(FASTA)
    if a.limit:
        ids, seqs = ids[:a.limit], seqs[:a.limit]
    print(f'{len(ids)} sequences, max len {max(len(s) for s in seqs)}',
          flush=True)

    flat = [(i, piece, n) for i, s in enumerate(seqs) for piece, n in chunks(s)]
    flat.sort(key=lambda t: len(t[1]))
    print(f'{len(flat)} chunks', flush=True)

    acc = np.zeros((len(ids), DIM), dtype=np.float64)
    wsum = np.zeros(len(ids), dtype=np.float64)

    checked = [False]

    def flush(batch):
        if not batch:
            return
        # batch_encode_plus was removed in transformers v5; __call__ is the
        # portable form and behaves identically.
        enc = tok([spaced(p) for _, p, _ in batch],
                  add_special_tokens=True, padding=True, return_tensors='pt')
        if not checked[0]:
            checked[0] = True
            n_tok = int(enc['attention_mask'][0].sum())
            n_res = len(batch[0][1])
            # ProtT5 appends </s> and prepends nothing, so residue j sits at
            # position j.  If a start token were added the offset would be
            # wrong by one and every embedding would be quietly misaligned.
            assert n_tok == n_res + 1, (
                f'expected {n_res}+1 tokens, got {n_tok}: the tokenizer is '
                f'adding a prefix token and rep[k, :len(p)] is misaligned')
            print(f'  token alignment OK ({n_res} residues -> {n_tok} tokens)',
                  flush=True)
        with torch.no_grad():
            rep = model(input_ids=enc['input_ids'].to(dev),
                        attention_mask=enc['attention_mask'].to(dev)
                        ).last_hidden_state
        for k, (i, p, n) in enumerate(batch):
            # no start token: residue j is at position j; </s> sits at len(p)
            v = rep[k, :len(p)].mean(0).float().cpu().numpy()
            acc[i] += v * n
            wsum[i] += n

    batch, ntok, done = [], 0, 0
    for item in flat:
        L = len(item[1]) + 1                # </s> only
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
    assert feats.shape[1] == DIM, f'expected {DIM} dims, got {feats.shape[1]}'
    assert np.isfinite(feats).all(), 'non-finite values in embeddings'

    os.makedirs(os.path.dirname(out_npy), exist_ok=True)
    np.save(out_npy, feats)
    json.dump(ids, open(out_ids, 'w'))
    print(f'wrote {feats.shape} to {out_npy}', flush=True)

    rng = np.random.default_rng(0)
    n = min(5000, len(ids))
    a_, b_ = rng.integers(0, len(ids), n), rng.integers(0, len(ids), n)
    cs = (feats[a_] * feats[b_]).sum(1) / (
        np.linalg.norm(feats[a_], axis=1) * np.linalg.norm(feats[b_], axis=1))
    print(f'mean cosine of random pairs = {cs.mean():.3f} '
          f'(collapse would be ~1)')

    # the ESM run must cover the same proteins, or the two arms are not
    # comparable and train_pdisjoint will silently zero-fill the difference
    esm_ids = f"{BASE}/notebooks/protein_ids_esm.json"
    if os.path.exists(esm_ids) and not a.limit:
        e = set(json.load(open(esm_ids)))
        p = set(ids)
        print(f'vs ESM ids: shared={len(e & p)} only_esm={len(e - p)} '
              f'only_prott5={len(p - e)}')


if __name__ == '__main__':
    main()
