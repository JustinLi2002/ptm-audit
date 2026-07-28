#!/usr/bin/env python3
"""
train_alldata.py — baseline / +PPI / permuted-PPI on the authors' final release.

The all_data.tar.gz release differs from the individually posted .tsv files in
three ways that matter here: methylation K and R are separate tasks, the test
partitions carry the natural class ratio, and the split is protein-level.
What it retains is the negative-sampling threshold — 35-97% of its training
proteins are still label-homogeneous.

That combination is exactly the cell of our design in which a protein-constant
feature channel is trained on homogeneous data and then evaluated on unseen
proteins. This script measures the effect directly on the released files rather
than on our reconstruction of them.

Validation is carved from the training partition at the protein level, so no
protein appears in both. Windows are 61 residues here rather than 31; the
architecture is unchanged because the recurrent state is length-independent.

Usage:
    python train_alldata.py --ptm acetylation_k --cond ppi --val-seed 0
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
DATA = f"{BASE}/allbig/all_data"
OUTDIR = f"{BASE}/alldata_runs"
PPI_FEAT = f"{BASE}/notebooks/protein_features_ppi.npy"
PPI_IDS = f"{BASE}/notebooks/protein_ids_ppi.json"

MAX_EPOCHS, PATIENCE = 100, 10
BATCH_TRAIN, BATCH_TEST, LR = 64, 512, 1e-3
VAL_FRAC = 0.10

# released file prefixes
PREFIX = {
    "acetylation_k": "acet_k", "glycosylation_n": "gly_n",
    "methylation_k": "met_k", "methylation_r": "met_r",
    "phosphorylation_st": "phos_st", "phosphorylation_y": "phos_y",
    "sumoylation_k": "sumo_k", "ubiquitination_k": "ubi_k",
}
PTMS = list(PREFIX)

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY") + ["U", "O", "X"]
aa_to_idx = {a: i for i, a in enumerate(AA_LIST)}
UNK = aa_to_idx["X"]
EYE = np.eye(len(AA_LIST), dtype=np.float32)


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
        idx = [aa_to_idx.get(a, UNK) for a in self.seq[i]]
        xs = torch.from_numpy(EYE[idx].T.astype(np.float32))
        if self.vecs is None:
            return xs, torch.tensor(self.y[i])
        return xs, torch.from_numpy(self.vecs[i]), torch.tensor(self.y[i])


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


def iqr_average(a):
    if a.shape[0] < 4:
        return a.mean(axis=0)
    q1, q3 = np.percentile(a, 25, axis=0), np.percentile(a, 75, axis=0)
    lo, hi = q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1)
    out = np.zeros(a.shape[1])
    for i in range(a.shape[1]):
        c = a[:, i]
        m = (c >= lo[i]) & (c <= hi[i])
        out[i] = c[m].mean() if m.any() else c.mean()
    return out


def train_one(trn_loader, val_loader, seed, device, vec_mode, ppi_dim, y_train):
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = (DeepMVP_PPI(ppi_dim=ppi_dim) if vec_mode else DeepMVP_Single()).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=LR)
    pw = torch.tensor([(y_train == 0).sum() / max((y_train == 1).sum(), 1)],
                      dtype=torch.float32, device=device)
    crit = nn.BCEWithLogitsLoss(pos_weight=pw)
    y_val = val_loader.dataset.y
    best, best_state, wait, ep = 0.0, None, 0, 0
    for ep in range(1, MAX_EPOCHS + 1):
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
        print(f"    epoch {ep:3d} val={va:.4f} best={best:.4f}", flush=True)
    model.load_state_dict(best_state)
    return model, best, ep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ptm", required=True, choices=PTMS)
    ap.add_argument("--cond", required=True, choices=["baseline", "ppi", "shuffled"])
    ap.add_argument("--val-seed", type=int, default=0)
    ap.add_argument("--n-models", type=int, default=2)
    ap.add_argument("--test-filter", default="90", choices=["70", "80", "90"])
    ap.add_argument("--shuffle-seed", type=int, default=42)
    a = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    pre = PREFIX[a.ptm]
    tag = f"{a.ptm}__{a.cond}__val{a.val_seed}__f{a.test_filter}"
    print(f"[{tag}] device={device}", flush=True)

    full = pd.read_csv(f"{DATA}/{pre}_training.tsv", sep="\t")
    tst = pd.read_csv(f"{DATA}/{pre}_testing_{a.test_filter}.tsv", sep="\t")

    # protein-level validation split, so train and val share no protein
    prots = np.array(sorted(full["protein"].unique()))
    rng = np.random.RandomState(a.val_seed)
    rng.shuffle(prots)
    n_val = max(1, int(len(prots) * VAL_FRAC))
    val_prots = set(prots[:n_val])
    val = full[full["protein"].isin(val_prots)].reset_index(drop=True)
    trn = full[~full["protein"].isin(val_prots)].reset_index(drop=True)

    overlap = set(trn["protein"]) & set(tst["protein"])
    print(f"  train={len(trn)} val={len(val)} test={len(tst)} "
          f"pos {trn.y.mean():.3f}/{val.y.mean():.3f}/{tst.y.mean():.3f} "
          f"| train-test protein overlap: {len(overlap)}", flush=True)
    for nm, d in (("train", trn), ("val", val), ("test", tst)):
        if len(d) == 0 or d.y.nunique() < 2:
            sys.exit(f"[{tag}] {nm} unusable (n={len(d)})")

    vec_mode = a.cond in ("ppi", "shuffled")
    ppi_dim = 128
    if vec_mode:
        feats = np.load(PPI_FEAT)
        with open(PPI_IDS) as fh:
            ids = json.load(fh)
        mapping = dict(zip(ids, feats))
        ppi_dim = feats.shape[1]
        if a.cond == "shuffled":
            r2 = np.random.RandomState(a.shuffle_seed)
            ps = list(mapping.keys())
            vs = [mapping[q] for q in ps]
            mapping = dict(zip(ps, [vs[i] for i in r2.permutation(len(ps))]))
            print(f"  permuted {len(ps)} vectors", flush=True)
        zero = np.zeros(ppi_dim, dtype=np.float32)
        cov = np.mean([p in mapping for p in tst["protein"]])
        print(f"  PPI coverage on test proteins: {100*cov:.1f}%", flush=True)
        vt, vv, vs_ = (np.stack([mapping.get(p, zero) for p in d["protein"]])
                       for d in (trn, val, tst))
    else:
        vt = vv = vs_ = None

    pin = device.type == "cuda"
    trn_loader = DataLoader(DS(trn, vt), BATCH_TRAIN, shuffle=True,
                            num_workers=4, pin_memory=pin, drop_last=True)
    val_loader = DataLoader(DS(val, vv), BATCH_TEST, shuffle=False, num_workers=2)
    tst_loader = DataLoader(DS(tst, vs_), BATCH_TEST, shuffle=False, num_workers=2)

    y_train, y_test = trn.y.values, tst.y.values
    probs, per_seed, epochs = [], [], []
    for s in range(a.n_models):
        print(f"  --- model seed {s}", flush=True)
        model, va, ep = train_one(trn_loader, val_loader, s, device,
                                  vec_mode, ppi_dim, y_train)
        pr = predict(model, tst_loader, device, vec_mode)
        probs.append(pr)
        per_seed.append(float(roc_auc_score(y_test, pr)))
        epochs.append(ep)
        print(f"  seed {s}: val={va:.4f} test={per_seed[-1]:.4f} (ep {ep})", flush=True)

    ens = iqr_average(np.stack(probs))
    res = dict(ptm=a.ptm, cond=a.cond, val_seed=a.val_seed,
               test_filter=a.test_filter, dataset="all_data",
               n_train=len(trn), n_val=len(val), n_test=len(tst),
               pos_train=float(trn.y.mean()), pos_test=float(tst.y.mean()),
               train_test_protein_overlap=len(overlap),
               auroc=float(roc_auc_score(y_test, ens)),
               auprc=float(average_precision_score(y_test, ens)),
               seed_aurocs=per_seed, seed_mean=float(np.mean(per_seed)),
               seed_sd=float(np.std(per_seed)), epochs=epochs)
    os.makedirs(OUTDIR, exist_ok=True)
    with open(f"{OUTDIR}/{tag}.json", "w") as fh:
        json.dump(res, fh, indent=2)
    out = tst[["protein", "aa", "pos", "y"]].copy()
    out["y_pred"] = ens
    out.to_csv(f"{OUTDIR}/{tag}.pred.tsv", sep="\t", index=False)
    print(f"[{tag}] AUROC={res['auroc']:.4f} AUPRC={res['auprc']:.4f} "
          f"seeds={res['seed_mean']:.4f}\u00b1{res['seed_sd']:.4f}", flush=True)


if __name__ == "__main__":
    main()
