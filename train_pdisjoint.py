#!/usr/bin/env python3
"""
train_pdisjoint.py — retrain under a protein-disjoint split.

Answers: does the *construction* of the training data affect generalisation to
unseen proteins? Two training sets differing only in the negative-sampling rule
(replica = >=10-site threshold, rebuilt = no threshold) are trained under the
same protein-disjoint split and evaluated on the SAME held-out test set, so the
only variable is how the training negatives were drawn.

Primary evaluation set: rebuilt_test (uniformly sampled negatives).
Secondary: replica_test, reported for reference only — its own construction
inflates scores, which is the point of the paper, not a result here.

One invocation = one (ptm, dataset, condition, split_seed) cell, training
--n-models models. Designed to be driven by a SLURM job array.

Usage:
    python train_pdisjoint.py --ptm methylation_r --dataset rebuilt \
        --cond baseline --split-seed 0 --n-models 2
"""

import argparse
import copy
import json
import os
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score, average_precision_score
from torch.utils.data import Dataset, DataLoader

BASE = "/home/FCAM/juli/HRP"
SPLIT_DIR = f"{BASE}/pdisjoint"
OUTDIR = f"{BASE}/pdisjoint_runs"
PPI_FEAT = f"{BASE}/notebooks/protein_features_ppi.npy"
PPI_IDS = f"{BASE}/notebooks/protein_ids_ppi.json"

MAX_EPOCHS, PATIENCE = 100, 10
BATCH_TRAIN, BATCH_TEST, LR = 64, 512, 1e-3

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY") + ["U", "O", "X"]
aa_to_idx = {a: i for i, a in enumerate(AA_LIST)}
EYE = np.eye(len(AA_LIST), dtype=np.float32)


# ── models: identical to the April/May training scripts ─────────────────────
class DeepMVP_Single(nn.Module):
    def __init__(self, in_ch=23, dropout=0.3):
        super().__init__()
        self.conv1 = nn.Conv1d(in_ch, 512, 5, padding=2); self.bn1 = nn.BatchNorm1d(512)
        self.conv2 = nn.Conv1d(512, 512, 5, padding=2);   self.bn2 = nn.BatchNorm1d(512)
        self.conv3 = nn.Conv1d(512, 128, 5, padding=2);   self.bn3 = nn.BatchNorm1d(128)
        self.drop = nn.Dropout(dropout)
        self.gru = nn.GRU(128, 50, batch_first=True, bidirectional=True)
        self.fc1 = nn.Linear(100, 64); self.bn4 = nn.BatchNorm1d(64)
        self.drop2 = nn.Dropout(dropout); self.fc_out = nn.Linear(64, 1)

    def forward(self, x):
        x = self.drop(F.leaky_relu(self.bn1(self.conv1(x))))
        x = self.drop(F.leaky_relu(self.bn2(self.conv2(x))))
        x = self.drop(F.leaky_relu(self.bn3(self.conv3(x))))
        _, h = self.gru(x.permute(0, 2, 1))
        x = torch.cat([h[0], h[1]], dim=1)
        x = self.drop2(F.leaky_relu(self.bn4(self.fc1(x))))
        return self.fc_out(x).squeeze(-1)


class DeepMVP_PPI(nn.Module):
    def __init__(self, in_ch=23, ppi_dim=128, dropout=0.3):
        super().__init__()
        self.conv1 = nn.Conv1d(in_ch, 512, 5, padding=2); self.bn1 = nn.BatchNorm1d(512)
        self.conv2 = nn.Conv1d(512, 512, 5, padding=2);   self.bn2 = nn.BatchNorm1d(512)
        self.conv3 = nn.Conv1d(512, 128, 5, padding=2);   self.bn3 = nn.BatchNorm1d(128)
        self.drop = nn.Dropout(dropout)
        self.gru = nn.GRU(128, 50, batch_first=True, bidirectional=True)
        self.ppi_fc = nn.Linear(ppi_dim, 128)
        self.fc1 = nn.Linear(228, 64); self.bn4 = nn.BatchNorm1d(64)
        self.drop2 = nn.Dropout(dropout); self.fc_out = nn.Linear(64, 1)

    def forward(self, x_seq, x_vec):
        x = self.drop(F.leaky_relu(self.bn1(self.conv1(x_seq))))
        x = self.drop(F.leaky_relu(self.bn2(self.conv2(x))))
        x = self.drop(F.leaky_relu(self.bn3(self.conv3(x))))
        _, h = self.gru(x.permute(0, 2, 1))
        x = torch.cat([h[0], h[1]], dim=1)
        p = F.leaky_relu(self.ppi_fc(x_vec))
        x = self.drop2(F.leaky_relu(self.bn4(self.fc1(torch.cat([x, p], 1)))))
        return self.fc_out(x).squeeze(-1)


class DS(Dataset):
    def __init__(self, df, vecs=None):
        self.seq = df["x"].values
        self.y = df["y"].values.astype(np.float32)
        self.vecs = vecs

    def __len__(self):
        return len(self.y)

    def __getitem__(self, i):
        idx = [aa_to_idx[a] for a in self.seq[i]]
        xs = torch.from_numpy(EYE[idx].T.astype(np.float32))
        if self.vecs is None:
            return xs, torch.tensor(self.y[i])
        return xs, torch.from_numpy(self.vecs[i]), torch.tensor(self.y[i])


def iqr_average(all_probs):
    if all_probs.shape[0] < 4:          # IQR trimming is meaningless below 4
        return all_probs.mean(axis=0)
    q1 = np.percentile(all_probs, 25, axis=0)
    q3 = np.percentile(all_probs, 75, axis=0)
    lo, hi = q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1)
    out = np.zeros(all_probs.shape[1])
    for i in range(all_probs.shape[1]):
        col = all_probs[:, i]
        m = (col >= lo[i]) & (col <= hi[i])
        out[i] = col[m].mean() if m.any() else col.mean()
    return out


def attach_vectors(df, mapping, dim):
    zero = np.zeros(dim, dtype=np.float32)
    return np.stack([mapping.get(p, zero) for p in df["protein"]])


def predict(model, loader, device, vec_mode):
    model.eval()
    out = []
    with torch.no_grad():
        for batch in loader:
            if vec_mode:
                xs, xv, _ = batch
                p = torch.sigmoid(model(xs.to(device), xv.to(device)))
            else:
                xs, _ = batch
                p = torch.sigmoid(model(xs.to(device)))
            out.append(p.cpu().numpy())
    return np.concatenate(out)


def train_one(trn_loader, val_loader, seed, device, vec_mode, ppi_dim, y_train):
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = (DeepMVP_PPI(ppi_dim=ppi_dim) if vec_mode else DeepMVP_Single()).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=LR)
    pw = torch.tensor([(y_train == 0).sum() / max((y_train == 1).sum(), 1)],
                      dtype=torch.float32, device=device)
    crit = nn.BCEWithLogitsLoss(pos_weight=pw)

    y_val = val_loader.dataset.y
    best, best_state, wait = 0.0, None, 0
    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        for batch in trn_loader:
            if vec_mode:
                xs, xv, yb = batch
                logits = model(xs.to(device), xv.to(device))
            else:
                xs, yb = batch
                logits = model(xs.to(device))
            opt.zero_grad(set_to_none=True)
            crit(logits, yb.to(device)).backward()
            opt.step()
        va = roc_auc_score(y_val, predict(model, val_loader, device, vec_mode))
        if va > best:
            best, best_state, wait = va, copy.deepcopy(model.state_dict()), 0
        else:
            wait += 1
            if wait >= PATIENCE:
                break
        print(f"    epoch {epoch:3d} val={va:.4f} best={best:.4f}", flush=True)
    model.load_state_dict(best_state)
    return model, best, epoch


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ptm", required=True)
    ap.add_argument("--dataset", required=True, choices=["replica", "rebuilt"])
    ap.add_argument("--cond", required=True,
                    choices=["baseline", "ppi", "shuffled"])
    ap.add_argument("--shuffle-seed", type=int, default=42,
                    help="seed for the protein-to-vector permutation")
    ap.add_argument("--split-seed", type=int, required=True)
    ap.add_argument("--n-models", type=int, default=2)
    ap.add_argument("--max-train", type=int, default=0,
                    help="subsample training rows (0 = use all)")
    a = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tag = f"{a.ptm}__{a.dataset}__{a.cond}__split{a.split_seed}"
    print(f"[{tag}] device={device}", flush=True)

    split = pd.read_csv(f"{SPLIT_DIR}/split_seed{a.split_seed}.csv")
    sp = dict(zip(split["accession"], split["split"]))

    train_src = pd.read_csv(f"{BASE}/{a.dataset}/{a.ptm}_all.tsv", sep="\t")
    train_src["split"] = train_src["protein"].map(sp)
    trn = train_src[train_src.split == "train"].reset_index(drop=True)
    val = train_src[train_src.split == "val"].reset_index(drop=True)

    # evaluation is always on the unrestricted-sampling test set
    eval_src = pd.read_csv(f"{BASE}/rebuilt/{a.ptm}_all.tsv", sep="\t")
    eval_src["split"] = eval_src["protein"].map(sp)
    tst = eval_src[eval_src.split == "test"].reset_index(drop=True)

    if a.max_train and len(trn) > a.max_train:
        trn = trn.sample(a.max_train, random_state=0).reset_index(drop=True)
        print(f"  subsampled train to {len(trn)}", flush=True)

    for nm, d in (("train", trn), ("val", val), ("test", tst)):
        if len(d) == 0 or d.y.nunique() < 2:
            sys.exit(f"[{tag}] {nm} split unusable (n={len(d)})")
    print(f"  train={len(trn)} val={len(val)} test={len(tst)} "
          f"pos: {trn.y.mean():.3f}/{val.y.mean():.3f}/{tst.y.mean():.3f}", flush=True)

    vec_mode = a.cond in ("ppi", "shuffled")
    ppi_dim = 128
    if vec_mode:
        feats = np.load(PPI_FEAT)
        with open(PPI_IDS) as fh:
            ids = json.load(fh)
        mapping = dict(zip(ids, feats))
        ppi_dim = feats.shape[1]
        if a.cond == "shuffled":
            # Permute which protein carries which vector. Preserves the marginal
            # distribution of vectors, within-protein constancy, and the set of
            # proteins that have no vector at all; destroys only the
            # protein-to-embedding correspondence. Under a protein-disjoint
            # split this asks whether the gain comes from structure in PPI space
            # or from something that survives randomisation.
            rng = np.random.RandomState(a.shuffle_seed)
            prots = list(mapping.keys())
            vecs = [mapping[q] for q in prots]
            order = rng.permutation(len(prots))
            mapping = dict(zip(prots, [vecs[i] for i in order]))
            n_moved = sum(1 for i, j in enumerate(order) if i != j)
            print(f"  permuted {len(prots)} vectors ({n_moved} moved)", flush=True)
        vt, vv, vs = (attach_vectors(d, mapping, ppi_dim) for d in (trn, val, tst))
    else:
        vt = vv = vs = None

    pin = device.type == "cuda"
    trn_loader = DataLoader(DS(trn, vt), BATCH_TRAIN, shuffle=True,
                            num_workers=4, pin_memory=pin, drop_last=True)
    val_loader = DataLoader(DS(val, vv), BATCH_TEST, shuffle=False, num_workers=2)
    tst_loader = DataLoader(DS(tst, vs), BATCH_TEST, shuffle=False, num_workers=2)

    y_train = trn.y.values
    y_test = tst.y.values
    probs, per_seed, epochs = [], [], []
    for s in range(a.n_models):
        print(f"  --- model seed {s}", flush=True)
        model, va, ep = train_one(trn_loader, val_loader, s, device,
                                  vec_mode, ppi_dim, y_train)
        pr = predict(model, tst_loader, device, vec_mode)
        probs.append(pr)
        per_seed.append(float(roc_auc_score(y_test, pr)))
        epochs.append(ep)
        print(f"  seed {s}: val={va:.4f} test={per_seed[-1]:.4f} (stopped ep {ep})",
              flush=True)

    ens = iqr_average(np.stack(probs))
    res = dict(
        ptm=a.ptm, dataset=a.dataset, cond=a.cond, split_seed=a.split_seed,
        n_train=len(trn), n_val=len(val), n_test=len(tst),
        pos_train=float(trn.y.mean()), pos_test=float(tst.y.mean()),
        auroc=float(roc_auc_score(y_test, ens)),
        auprc=float(average_precision_score(y_test, ens)),
        seed_aurocs=per_seed, seed_mean=float(np.mean(per_seed)),
        seed_sd=float(np.std(per_seed)), epochs=epochs,
    )
    os.makedirs(OUTDIR, exist_ok=True)
    with open(f"{OUTDIR}/{tag}.json", "w") as fh:
        json.dump(res, fh, indent=2)
    out = tst[["protein", "aa", "pos", "y"]].copy()
    out["y_pred"] = ens
    out.to_csv(f"{OUTDIR}/{tag}.pred.tsv", sep="\t", index=False)
    print(f"[{tag}] AUROC={res['auroc']:.4f} AUPRC={res['auprc']:.4f} "
          f"seeds={res['seed_mean']:.4f}±{res['seed_sd']:.4f}", flush=True)


if __name__ == "__main__":
    main()
