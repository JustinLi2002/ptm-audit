#!/usr/bin/env python3
"""
rerun_inference.py — regenerate per-site test predictions from saved checkpoints.

Loads the checkpoints trained in April-May 2026 and re-runs forward passes on the
test sets, producing:
  * per-site y_pred TSVs (needed for the restricted-subset comparison)
  * per-seed test AUROC  (needed for error bars — the logs only recorded val AUC)
  * IQR-ensemble AUROC   (self-check against Table 2 of the manuscript)
  * restricted-subset AUROC (label-homogeneous proteins removed)

No retraining. Model definitions are copied verbatim from
ContextPTM/models/deepmvp_reproduce_v2.py and deepmvp_ppi.py; the script asserts
that every checkpoint's state_dict matches before running, so any divergence
fails loudly instead of silently producing wrong numbers.

Usage:
    python rerun_inference.py --cond baseline
    python rerun_inference.py --cond ppi --ptm methylation_k
    python rerun_inference.py --cond all
"""

import argparse
import collections
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

BASE = os.environ.get("PTM_AUDIT_BASE", "/home/FCAM/juli/HRP")
RETRAIN = f"{BASE}/retrain"
OUTDIR = f"{BASE}/inference_out"
BATCH_TEST = 512

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY") + ["U", "O", "X"]
aa_to_idx = {aa: i for i, aa in enumerate(AA_LIST)}
EYE = np.eye(len(AA_LIST), dtype=np.float32)

PTMS = ["acetylation_k", "glycosylation_n", "methylation_k", "methylation_r",
        "phosphorylation_st", "phosphorylation_y", "sumoylation_k",
        "ubiquitination_k"]

# condition -> (checkpoint dir, feature .npy, ids .json)
CONDS = {
    "baseline": (f"{BASE}/checkpoints", None, None),
    "kinase":   (f"{BASE}/checkpoints_kinase",
                 f"{BASE}/notebooks/protein_features.npy",
                 f"{BASE}/notebooks/protein_ids.json"),
    "ppi":      (f"{BASE}/checkpoints_ppi",
                 f"{BASE}/notebooks/protein_features_ppi.npy",
                 f"{BASE}/notebooks/protein_ids_ppi.json"),
    "shuffled": (f"{BASE}/checkpoints_kinase_shuffled/shuffled",
                 f"{BASE}/notebooks/protein_features.npy",
                 f"{BASE}/notebooks/protein_ids.json"),
}

# manuscript Table 2, for the self-check: (baseline, kinase, ppi)
TABLE2 = {
    "acetylation_k":      (0.9049, 0.9615, 0.9560),
    "glycosylation_n":    (0.9868, 0.9968, 0.9954),
    "methylation_k":      (0.9503, 0.9814, 0.9785),
    "methylation_r":      (0.9306, 0.9810, 0.9738),
    "phosphorylation_st": (0.9510, None,   0.9654),
    "phosphorylation_y":  (0.8710, None,   0.9220),
    "sumoylation_k":      (0.8615, 0.9256, 0.9210),
    "ubiquitination_k":   (0.8814, 0.9296, 0.9248),
}


# ── models (verbatim from the training scripts) ──────────────────────────────
class DeepMVP_Single(nn.Module):
    def __init__(self, in_ch=23, seq_len=31, dropout=0.3):
        super().__init__()
        self.conv1 = nn.Conv1d(in_ch, 512, kernel_size=5, padding=2)
        self.bn1 = nn.BatchNorm1d(512)
        self.conv2 = nn.Conv1d(512, 512, kernel_size=5, padding=2)
        self.bn2 = nn.BatchNorm1d(512)
        self.conv3 = nn.Conv1d(512, 128, kernel_size=5, padding=2)
        self.bn3 = nn.BatchNorm1d(128)
        self.drop = nn.Dropout(dropout)
        self.gru = nn.GRU(128, 50, batch_first=True, bidirectional=True)
        self.fc1 = nn.Linear(100, 64)
        self.bn4 = nn.BatchNorm1d(64)
        self.drop2 = nn.Dropout(dropout)
        self.fc_out = nn.Linear(64, 1)

    def forward(self, x):
        x = self.drop(F.leaky_relu(self.bn1(self.conv1(x))))
        x = self.drop(F.leaky_relu(self.bn2(self.conv2(x))))
        x = self.drop(F.leaky_relu(self.bn3(self.conv3(x))))
        x = x.permute(0, 2, 1)
        _, h = self.gru(x)
        x = torch.cat([h[0], h[1]], dim=1)
        x = self.drop2(F.leaky_relu(self.bn4(self.fc1(x))))
        return self.fc_out(x).squeeze(-1)


class DeepMVP_PPI(nn.Module):
    def __init__(self, in_ch=23, ppi_dim=605, dropout=0.3):
        super().__init__()
        self.conv1 = nn.Conv1d(in_ch, 512, kernel_size=5, padding=2)
        self.bn1 = nn.BatchNorm1d(512)
        self.conv2 = nn.Conv1d(512, 512, kernel_size=5, padding=2)
        self.bn2 = nn.BatchNorm1d(512)
        self.conv3 = nn.Conv1d(512, 128, kernel_size=5, padding=2)
        self.bn3 = nn.BatchNorm1d(128)
        self.drop = nn.Dropout(dropout)
        self.gru = nn.GRU(128, 50, batch_first=True, bidirectional=True)
        self.ppi_fc = nn.Linear(ppi_dim, 128)
        self.fc1 = nn.Linear(228, 64)
        self.bn4 = nn.BatchNorm1d(64)
        self.drop2 = nn.Dropout(dropout)
        self.fc_out = nn.Linear(64, 1)

    def forward(self, x_seq, x_vec):
        x = self.drop(F.leaky_relu(self.bn1(self.conv1(x_seq))))
        x = self.drop(F.leaky_relu(self.bn2(self.conv2(x))))
        x = self.drop(F.leaky_relu(self.bn3(self.conv3(x))))
        x = x.permute(0, 2, 1)
        _, h = self.gru(x)
        x = torch.cat([h[0], h[1]], dim=1)
        p = F.leaky_relu(self.ppi_fc(x_vec))
        x = torch.cat([x, p], dim=1)
        x = self.drop2(F.leaky_relu(self.bn4(self.fc1(x))))
        return self.fc_out(x).squeeze(-1)


class DS_Seq(Dataset):
    def __init__(self, df):
        self.seq = df["x"].values
        self.y = df["y"].values.astype(np.float32)

    def __len__(self):
        return len(self.y)

    def __getitem__(self, i):
        idx = [aa_to_idx[a] for a in self.seq[i]]
        return torch.from_numpy(EYE[idx].T.astype(np.float32)), torch.tensor(self.y[i])


class DS_Vec(Dataset):
    def __init__(self, df):
        self.seq = df["x"].values
        self.vec = df["VEC"].values
        self.y = df["y"].values.astype(np.float32)

    def __len__(self):
        return len(self.y)

    def __getitem__(self, i):
        idx = [aa_to_idx[a] for a in self.seq[i]]
        return (torch.from_numpy(EYE[idx].T.astype(np.float32)),
                torch.from_numpy(np.asarray(self.vec[i], dtype=np.float32)),
                torch.tensor(self.y[i]))


def iqr_average(all_probs):
    q1 = np.percentile(all_probs, 25, axis=0)
    q3 = np.percentile(all_probs, 75, axis=0)
    iqr = q3 - q1
    lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    out = np.zeros(all_probs.shape[1])
    for i in range(all_probs.shape[1]):
        col = all_probs[:, i]
        m = (col >= lo[i]) & (col <= hi[i])
        out[i] = col[m].mean() if m.any() else col.mean()
    return out


def load_features(feat_path, ids_path, seed_shuffle=None):
    feats = np.load(feat_path)
    with open(ids_path) as fh:
        ids = json.load(fh)
    mapping = dict(zip(ids, feats))
    if seed_shuffle is not None:
        # reproduce make_shuffled_mapping() exactly
        rng = np.random.RandomState(seed_shuffle)
        prots = list(mapping.keys())
        vecs = [mapping[p] for p in prots]
        mapping = dict(zip(prots, [vecs[i] for i in rng.permutation(len(prots))]))
    return mapping, feats.shape[1]


def homogeneous_proteins(ptm):
    per = collections.defaultdict(lambda: [0, 0])
    df = pd.read_csv(f"{RETRAIN}/{ptm}_train.tsv", sep="\t")
    for p, y in zip(df["protein"], df["y"]):
        per[p][int(y)] += 1
    return {p for p, (neg, pos) in per.items() if neg == 0 or pos == 0}


def run(cond, ptm, device):
    ckpt_dir, feat_path, ids_path = CONDS[cond]
    d = f"{ckpt_dir}/{ptm}"
    if not os.path.isdir(d):
        print(f"  [{cond}/{ptm}] no checkpoint dir, skipped")
        return None

    test = pd.read_csv(f"{RETRAIN}/{ptm}_test.tsv", sep="\t")

    if feat_path is None:
        model_fn, ds = (lambda: DeepMVP_Single()), DS_Seq(test)
        vec_mode = False
    else:
        mapping, dim = load_features(
            feat_path, ids_path, seed_shuffle=42 if cond == "shuffled" else None)
        zero = np.zeros(dim, dtype=np.float32)
        test = test.copy()
        test["VEC"] = [mapping.get(p, zero) for p in test["protein"]]
        model_fn, ds = (lambda: DeepMVP_PPI(ppi_dim=dim)), DS_Vec(test)
        vec_mode = True

    loader = DataLoader(ds, BATCH_TEST, shuffle=False, num_workers=2)
    y_true = test["y"].values.astype(int)

    probs, per_seed = [], []
    for seed in range(10):
        f = f"{d}/model_{seed}.pt"
        if not os.path.exists(f):
            print(f"  [{cond}/{ptm}] missing model_{seed}.pt")
            continue
        model = model_fn().to(device)
        state = torch.load(f, map_location=device)
        missing, unexpected = model.load_state_dict(state, strict=False)
        if missing or unexpected:
            sys.exit(f"STATE DICT MISMATCH {f}\n  missing={missing}\n  "
                     f"unexpected={unexpected}\n  architecture has diverged — stop.")
        model.eval()
        ps = []
        with torch.no_grad():
            for batch in loader:
                if vec_mode:
                    xs, xv, _ = batch
                    p = torch.sigmoid(model(xs.to(device), xv.to(device)))
                else:
                    xs, _ = batch
                    p = torch.sigmoid(model(xs.to(device)))
                ps.append(p.cpu().numpy())
        pr = np.concatenate(ps)
        probs.append(pr)
        per_seed.append(roc_auc_score(y_true, pr))

    if not probs:
        return None

    ens = iqr_average(np.stack(probs))
    auc = roc_auc_score(y_true, ens)
    aup = average_precision_score(y_true, ens)

    homo = homogeneous_proteins(ptm)
    mask = ~test["protein"].isin(homo).values
    auc_r = roc_auc_score(y_true[mask], ens[mask]) if mask.sum() and \
        len(set(y_true[mask])) > 1 else float("nan")

    os.makedirs(OUTDIR, exist_ok=True)
    out = test[["protein", "aa", "pos", "y"]].copy()
    out["y_pred"] = ens
    for i, pr in enumerate(probs):
        out[f"seed_{i}"] = pr
    out.to_csv(f"{OUTDIR}/{cond}__{ptm}.tsv", sep="\t", index=False)

    ref = TABLE2[ptm][{"baseline": 0, "kinase": 1, "ppi": 2,
                       "shuffled": 1}[cond]]
    flag = ""
    if ref is not None:
        flag = " OK" if abs(auc - ref) < 0.005 else f" MISMATCH (Table2={ref:.4f})"
    print(f"  {cond:9s} {ptm:20s} ens={auc:.4f} auprc={aup:.4f} "
          f"seeds={np.mean(per_seed):.4f}±{np.std(per_seed):.4f} "
          f"restricted={auc_r:.4f} n_r={mask.sum()}{flag}")
    return dict(cond=cond, ptm=ptm, ensemble=auc, auprc=aup,
                seed_mean=float(np.mean(per_seed)), seed_sd=float(np.std(per_seed)),
                seed_aurocs=[float(x) for x in per_seed],
                restricted=float(auc_r), n_restricted=int(mask.sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cond", default="all",
                    choices=list(CONDS) + ["all"])
    ap.add_argument("--ptm", default="all")
    a = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device: {device}")
    conds = list(CONDS) if a.cond == "all" else [a.cond]
    ptms = PTMS if a.ptm == "all" else [a.ptm]

    results = []
    for c in conds:
        for p in ptms:
            r = run(c, p, device)
            if r:
                results.append(r)

    os.makedirs(OUTDIR, exist_ok=True)
    with open(f"{OUTDIR}/summary.json", "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nwrote {len(results)} results to {OUTDIR}/summary.json")
    bad = [r for r in results if TABLE2[r["ptm"]][
        {"baseline": 0, "kinase": 1, "ppi": 2, "shuffled": 1}[r["cond"]]] is not None
        and abs(r["ensemble"] - TABLE2[r["ptm"]][
            {"baseline": 0, "kinase": 1, "ppi": 2, "shuffled": 1}[r["cond"]]]) >= 0.005]
    if bad:
        print(f"\n!! {len(bad)} results disagree with Table 2 by >=0.005 — "
              f"investigate before using any of this.")


if __name__ == "__main__":
    main()
