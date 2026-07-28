#!/usr/bin/env python3
"""
Protein-disjoint train/val/test split for ContextPTM.

Input:
  --clstr   CD-HIT .clstr file (run at -c 0.4 -n 2 over all proteins in the dataset)
  --sites   CSV/TSV with at least: accession, ptm_type, label   (one row per candidate site)
Output:
  --out     CSV: accession, cluster_id, split

Key property: ONE global split shared by all PTM types. A protein assigned to
test for Phospho S/T is in test for every other PTM type too. Without this the
cross-PTM transfer matrix is still leaky.
"""

import argparse
import random
import re
import sys
from collections import defaultdict

import pandas as pd

# matches ">sp|P12345|NAME_HUMAN" or ">P12345" inside a .clstr line
ACC_RE = re.compile(r">(?:\w{2}\|)?([A-Za-z0-9_.\-]+)")


def parse_clstr(path):
    """Return {accession: cluster_id}. Accessions are stripped of |NAME suffix."""
    mapping = {}
    cluster_id = None
    n_lines = 0
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">Cluster"):
                cluster_id = int(line.split()[1])
                continue
            m = ACC_RE.search(line)
            if m is None:
                print(f"WARN: unparsed .clstr line: {line[:80]}", file=sys.stderr)
                continue
            acc = m.group(1).rstrip('.')
            # sp|P12345|NAME_HUMAN -> regex already grabbed P12345 via \w{2}\| prefix
            if acc in mapping:
                raise ValueError(f"accession {acc} appears in two clusters")
            mapping[acc] = cluster_id
            n_lines += 1
    print(f"parsed {n_lines} sequences in {cluster_id + 1} clusters", file=sys.stderr)
    return mapping


def assign_splits(site_df, acc2cluster, test_frac, val_frac, seed):
    """Greedy: shuffle clusters, fill test then val by cumulative site count."""
    missing = set(site_df["accession"]) - set(acc2cluster)
    if missing:
        raise SystemExit(
            f"ERROR: {len(missing)} accessions in the site table are absent from the "
            f".clstr file (e.g. {sorted(missing)[:5]}). The FASTA fed to CD-HIT must "
            f"cover every protein in the dataset."
        )

    sites_per_cluster = defaultdict(int)
    for acc, n in site_df.groupby("accession").size().items():
        sites_per_cluster[acc2cluster[acc]] += n

    clusters = sorted(sites_per_cluster)  # sort first so seed fully determines order
    rng = random.Random(seed)
    rng.shuffle(clusters)

    total = sum(sites_per_cluster.values())
    test_target = total * test_frac
    val_target = total * val_frac

    split_of_cluster, acc_test, acc_val = {}, 0, 0
    for c in clusters:
        n = sites_per_cluster[c]
        if acc_test < test_target:
            split_of_cluster[c] = "test"
            acc_test += n
        elif acc_val < val_target:
            split_of_cluster[c] = "val"
            acc_val += n
        else:
            split_of_cluster[c] = "train"
    return split_of_cluster


def report(site_df, out_df):
    merged = site_df.merge(out_df, on="accession", how="left")

    # hard assertion: no protein straddles two splits
    straddle = merged.groupby("accession")["split"].nunique()
    assert (straddle == 1).all(), "protein assigned to multiple splits"

    print("\n=== sites per split, per PTM type ===")
    tab = merged.pivot_table(
        index="ptm_type", columns="split", values="label", aggfunc="size", fill_value=0
    )
    print(tab.div(tab.sum(axis=1), axis=0).round(3))

    print("\n=== positive rate per split, per PTM type ===")
    print(
        merged.pivot_table(
            index="ptm_type", columns="split", values="label", aggfunc="mean"
        ).round(4)
    )

    print("\n=== absolute test-set size (watch for tiny cells) ===")
    print(tab["test"] if "test" in tab else "no test column")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--clstr", required=True)
    p.add_argument("--sites", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--test-frac", type=float, default=0.10)
    p.add_argument("--val-frac", type=float, default=0.10)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--sep", default=",")
    a = p.parse_args()

    site_df = pd.read_csv(a.sites, sep=a.sep)
    for col in ("accession", "ptm_type", "label"):
        if col not in site_df.columns:
            raise SystemExit(f"missing required column: {col}")

    acc2cluster = parse_clstr(a.clstr)
    split_of_cluster = assign_splits(
        site_df, acc2cluster, a.test_frac, a.val_frac, a.seed
    )

    known = set(site_df["accession"])
    out_df = pd.DataFrame(
        [
            {"accession": acc, "cluster_id": c, "split": split_of_cluster[c]}
            for acc, c in acc2cluster.items()
            if acc in known
        ]
    )
    out_df.to_csv(a.out, index=False)
    report(site_df, out_df)
    print(f"\nwrote {a.out}  ({len(out_df)} proteins)")


if __name__ == "__main__":
    main()
