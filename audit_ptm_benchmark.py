#!/usr/bin/env python3
"""
audit_ptm_benchmark.py — label-leakage audit for multi-PTM site prediction benchmarks.

Input: a directory of TSVs named <ptm>_{train,test}.tsv with columns
       protein / aa / pos / x / y   (x = sequence window, y = 0/1 label)

Reports:
  T1  dataset size, class balance, test fraction
  T2  residue composition + centre-residue self-consistency
  T3  label-homogeneous proteins (the leakage channel)
  T4  lookup-rule accuracy on the test set
  T5  protein-identity-only baseline AUROC vs a sequence model
  T6  cross-dataset site overlap

Usage:
    python audit_ptm_benchmark.py --dir retrain
    python audit_ptm_benchmark.py --dir retrain --seq-auroc seq_auroc.json
"""

import argparse
import collections
import glob
import itertools
import json
import os
import sys


def read_split(path):
    """Yield (protein, aa, pos, window, label). Tolerates missing window column."""
    with open(path) as fh:
        header = next(fh).rstrip("\n").split("\t")
        try:
            iy = header.index("y")
        except ValueError:
            sys.exit(f"{path}: no 'y' column; found {header}")
        ix = header.index("x") if "x" in header else None
        for line in fh:
            p = line.rstrip("\n").split("\t")
            yield p[0], p[1], p[2], (p[ix] if ix is not None else ""), int(p[iy])


def discover(d):
    """Return sorted PTM names that have both a train and a test file."""
    names = set()
    for f in glob.glob(os.path.join(d, "*_train.tsv")):
        names.add(os.path.basename(f)[: -len("_train.tsv")])
    return sorted(n for n in names if os.path.exists(os.path.join(d, f"{n}_test.tsv")))


def auroc(pairs):
    """Rank-based AUROC with correct handling of ties (mid-ranks)."""
    pairs = sorted(pairs, key=lambda t: t[0])
    P = sum(l for _, l in pairs)
    N = len(pairs) - P
    if P == 0 or N == 0:
        return float("nan")
    rank_sum = 0.0
    i = 0
    while i < len(pairs):
        j = i
        while j < len(pairs) and pairs[j][0] == pairs[i][0]:
            j += 1
        mid = (i + j + 1) / 2.0
        rank_sum += sum(mid for k in range(i, j) if pairs[k][1] == 1)
        i = j
    return (rank_sum - P * (P + 1) / 2.0) / (P * N)


def load(d, names):
    """data[ptm][split] = list of records."""
    data = {}
    for nm in names:
        data[nm] = {
            sp: list(read_split(os.path.join(d, f"{nm}_{sp}.tsv")))
            for sp in ("train", "test")
        }
    return data


def t1_balance(data):
    print("\n=== T1  size / class balance / test fraction ===")
    print(f"{'PTM':22s} {'train':>8s} {'tr pos%':>8s} {'test':>7s} {'te pos%':>8s}"
          f" {'delta':>7s} {'test frac':>9s}")
    for nm, d in data.items():
        tr, te = d["train"], d["test"]
        trp = sum(r[4] for r in tr) / len(tr)
        tep = sum(r[4] for r in te) / len(te)
        print(f"{nm:22s} {len(tr):8d} {trp:8.3f} {len(te):7d} {tep:8.3f}"
              f" {tep - trp:+7.3f} {len(te) / (len(tr) + len(te)):9.3f}")
    print("  NOTE: a test positive rate pinned at exactly 0.500 indicates a")
    print("        deliberately balanced test set, not stratified sampling.")


def t2_residues(data):
    print("\n=== T2  residue composition / centre-residue consistency ===")
    print(f"{'PTM':22s} {'split':>6s} {'n':>8s} {'centre mismatch':>16s}  composition")
    for nm, d in data.items():
        for sp, recs in d.items():
            comp = collections.Counter(r[1] for r in recs)
            bad = sum(1 for r in recs
                      if r[3] and (len(r[3]) % 2 == 0 or r[3][len(r[3]) // 2] != r[1]))
            print(f"{nm:22s} {sp:>6s} {len(recs):8d} "
                  f"{bad:6d} ({100 * bad / len(recs):5.2f}%) "
                  f"  {dict(comp.most_common(4))}")
    print("  NOTE: a file whose composition disagrees with its name (e.g. a")
    print("        '*_k' file containing R) is a merged task, not a pure one.")


def homogeneous_map(train_recs):
    """{protein: 1|0} for proteins whose training labels are all-positive/all-negative."""
    per = collections.defaultdict(lambda: [0, 0])
    for prot, _, _, _, y in train_recs:
        per[prot][y] += 1
    return per, {p: (1 if neg == 0 else 0)
                 for p, (neg, pos) in per.items() if neg == 0 or pos == 0}


def t3_homogeneity(data):
    print("\n=== T3  label-homogeneous proteins in TRAIN (the leakage channel) ===")
    print(f"{'PTM':22s} {'proteins':>9s} {'pure+':>7s} {'pure-':>7s}"
          f" {'pure+ %prot':>12s} {'pure+ %sites':>13s} {'pure- %sites':>13s}")
    out = {}
    for nm, d in data.items():
        per, _ = homogeneous_map(d["train"])
        tot = sum(neg + pos for neg, pos in per.values())
        pp = [v for v in per.values() if v[0] == 0]
        pn = [v for v in per.values() if v[1] == 0]
        sp = sum(pos for neg, pos in pp)
        sn = sum(neg for neg, pos in pn)
        out[nm] = 100 * sp / tot
        print(f"{nm:22s} {len(per):9d} {len(pp):7d} {len(pn):7d}"
              f" {100 * len(pp) / len(per):11.1f}% {100 * sp / tot:12.1f}%"
              f" {100 * sn / tot:12.1f}%")
    print("  NOTE: pure-positive proteins are a CONSTRUCTIVE artefact when negatives")
    print("        are sampled only from proteins above a site-count threshold.")
    return out


def t4_lookup(data):
    print("\n=== T4  lookup rule on TEST (no sequence information used) ===")
    print(f"{'PTM':22s} {'test':>8s} {'lookupable':>11s} {'frac':>7s} {'accuracy':>9s}")
    for nm, d in data.items():
        _, homo = homogeneous_map(d["train"])
        n = hit = corr = 0
        for prot, _, _, _, y in d["test"]:
            n += 1
            if prot in homo:
                hit += 1
                corr += (homo[prot] == y)
        acc = 100 * corr / hit if hit else float("nan")
        print(f"{nm:22s} {n:8d} {hit:11d} {100 * hit / n:6.1f}% {acc:8.1f}%")
    print("  NOTE: accuracy at 100% confirms the channel is constructive rather")
    print("        than coincidental; ~80% indicates a chance-driven subset.")


def t5_identity_auroc(data, seq_auroc):
    print("\n=== T5  protein-identity-only baseline AUROC ===")
    print(f"{'PTM':22s} {'identity':>9s} {'sequence':>9s} {'ratio':>7s}")
    for nm, d in data.items():
        per, _ = homogeneous_map(d["train"])
        # the most a protein-constant feature can encode: per-protein positive rate
        score = {p: pos / (neg + pos) for p, (neg, pos) in per.items()}
        pairs = [(score.get(prot, 0.5), y) for prot, _, _, _, y in d["test"]]
        a = auroc(pairs)
        s = seq_auroc.get(nm)
        if s:
            print(f"{nm:22s} {a:9.4f} {s:9.4f} {a / s:7.2f}")
        else:
            print(f"{nm:22s} {a:9.4f} {'-':>9s} {'-':>7s}")
    print("  NOTE: this classifier reads zero amino acids. Any protein-constant")
    print("        feature channel can in principle recover it.")


def t6_overlap(data):
    print("\n=== T6  cross-dataset site overlap (train_i vs test_j) ===")
    sites = {nm: {sp: {(r[0], r[2]) for r in recs} for sp, recs in d.items()}
             for nm, d in data.items()}
    flagged = False
    for a, b in itertools.permutations(sites, 2):
        inter = sites[a]["train"] & sites[b]["test"]
        if inter:
            pct = 100 * len(inter) / len(sites[b]["test"])
            if pct >= 1.0:
                flagged = True
                print(f"  {a}_train  covers {pct:5.1f}% of {b}_test "
                      f"({len(inter)}/{len(sites[b]['test'])})")
    if not flagged:
        print("  no cross-dataset overlap above 1%")
    else:
        print("  NOTE: overlap here invalidates the corresponding cross-PTM")
        print("        transfer cells — those are memorisation, not generalisation.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="directory of *_{train,test}.tsv")
    ap.add_argument("--seq-auroc", help="JSON {ptm_name: auroc} for the T5 comparison")
    ap.add_argument("--skip", nargs="*", default=[], help="table ids to skip, e.g. T6")
    a = ap.parse_args()

    names = discover(a.dir)
    if not names:
        sys.exit(f"no *_train.tsv / *_test.tsv pairs found in {a.dir}")
    print(f"auditing {len(names)} PTM types from {a.dir}: {', '.join(names)}")

    seq = {}
    if a.seq_auroc:
        with open(a.seq_auroc) as fh:
            seq = json.load(fh)

    data = load(a.dir, names)
    if "T1" not in a.skip: t1_balance(data)
    if "T2" not in a.skip: t2_residues(data)
    if "T3" not in a.skip: t3_homogeneity(data)
    if "T4" not in a.skip: t4_lookup(data)
    if "T5" not in a.skip: t5_identity_auroc(data, seq)
    if "T6" not in a.skip: t6_overlap(data)


if __name__ == "__main__":
    main()
