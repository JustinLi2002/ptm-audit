#!/usr/bin/env python3
"""
restricted_eval.py — the decisive comparison for paper 1.

For every site_prediction.tsv under runs/, this script
  1. identifies which test set it belongs to BY CONTENT (never by directory name),
  2. computes AUROC on the full test set,
  3. recomputes it after removing sites on label-homogeneous training proteins,
so the sequence model and the identity baseline can be compared on identical
subsets.

Usage:
    python restricted_eval.py --runs deepmvp/DeepMVP/runs --data retrain
    python restricted_eval.py --runs ... --data ... --diagonal-only
"""

import argparse
import collections
import glob
import os
import sys

# Table 2 of the manuscript, for identifying which condition a run corresponds to
TABLE2 = {
    'phosphorylation_st': (0.9510, 0.9654), 'phosphorylation_y': (0.8710, 0.9220),
    'sumoylation_k':      (0.8615, 0.9210), 'ubiquitination_k':  (0.8814, 0.9248),
    'acetylation_k':      (0.9049, 0.9560), 'glycosylation_n':   (0.9868, 0.9954),
    'methylation_k':      (0.9503, 0.9785), 'methylation_r':     (0.9306, 0.9738),
}
# identity-baseline AUROCs already computed (full / homogeneous-removed)
IDENTITY = {
    'acetylation_k': (0.8765, 0.7494), 'glycosylation_n': (0.9247, 0.5914),
    'methylation_k': (0.9571, 0.7830), 'methylation_r':   (0.9450, 0.8001),
    'phosphorylation_st': (0.7729, 0.7300), 'phosphorylation_y': (0.8020, 0.7527),
    'sumoylation_k': (0.8201, 0.7092), 'ubiquitination_k': (0.7990, 0.7229),
}


def auroc(pairs):
    pairs = sorted(pairs, key=lambda t: t[0])
    P = sum(l for _, l in pairs)
    N = len(pairs) - P
    if P == 0 or N == 0:
        return float('nan')
    rs = 0.0
    i = 0
    while i < len(pairs):
        j = i
        while j < len(pairs) and pairs[j][0] == pairs[i][0]:
            j += 1
        mid = (i + j + 1) / 2.0
        rs += sum(mid for k in range(i, j) if pairs[k][1] == 1)
        i = j
    return (rs - P * (P + 1) / 2.0) / (P * N)


def read_cols(path, need):
    """Return list of dicts for the requested columns."""
    with open(path) as fh:
        head = next(fh).rstrip('\n').split('\t')
        idx = {}
        for c in need:
            if c not in head:
                return None
            idx[c] = head.index(c)
        rows = []
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) <= max(idx.values()):
                continue
            rows.append({c: p[i] for c, i in idx.items()})
        return rows


def build_reference(data_dir):
    """test-set site fingerprints + homogeneous-protein sets from train."""
    ref = {}
    for nm in TABLE2:
        te = read_cols(os.path.join(data_dir, f'{nm}_test.tsv'), ['protein', 'pos'])
        tr = read_cols(os.path.join(data_dir, f'{nm}_train.tsv'),
                       ['protein', 'pos', 'y'])
        if te is None or tr is None:
            sys.exit(f'missing or malformed {nm} files in {data_dir}')
        per = collections.defaultdict(lambda: [0, 0])
        for r in tr:
            per[r['protein']][int(r['y'])] += 1
        ref[nm] = {
            'sites': {(r['protein'], r['pos']) for r in te},
            'homo': {p for p, (neg, pos) in per.items() if neg == 0 or pos == 0},
        }
    return ref


def identify(sites, ref):
    """Return the PTM whose test set best matches this prediction file."""
    best, best_j = None, 0.0
    for nm, r in ref.items():
        inter = len(sites & r['sites'])
        union = len(sites | r['sites'])
        j = inter / union if union else 0.0
        if j > best_j:
            best, best_j = nm, j
    return best, best_j


def label_condition(nm, a):
    """Guess Baseline vs +PPI by proximity to the published Table 2 values."""
    b, p = TABLE2[nm]
    db, dp = abs(a - b), abs(a - p)
    if min(db, dp) > 0.02:
        return 'unmatched'
    return 'Baseline' if db < dp else '+PPI'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', required=True)
    ap.add_argument('--data', required=True)
    ap.add_argument('--diagonal-only', action='store_true',
                    help='only report files whose parent dir names agree (self-model)')
    a = ap.parse_args()

    ref = build_reference(a.data)
    files = sorted(glob.glob(os.path.join(a.runs, '**', 'site_prediction.tsv'),
                             recursive=True))
    if not files:
        sys.exit(f'no site_prediction.tsv under {a.runs}')
    print(f'found {len(files)} prediction files\n')

    hdr = (f"{'file':52s} {'testset':20s} {'match':>6s} {'n':>7s} "
           f"{'full':>7s} {'cond':>9s} {'restricted':>11s} {'id-restr':>9s}")
    print(hdr)
    print('-' * len(hdr))

    for f in files:
        rows = read_cols(f, ['protein', 'pos', 'y', 'y_pred'])
        if rows is None:
            print(f'{os.path.relpath(f, a.runs):52s} (no y_pred column, skipped)')
            continue
        sites = {(r['protein'], r['pos']) for r in rows}
        nm, jac = identify(sites, ref)
        if nm is None or jac < 0.5:
            print(f'{os.path.relpath(f, a.runs):52s} UNIDENTIFIED (best J={jac:.2f})')
            continue

        homo = ref[nm]['homo']
        full = [(float(r['y_pred']), int(r['y'])) for r in rows]
        rest = [(float(r['y_pred']), int(r['y'])) for r in rows
                if r['protein'] not in homo]
        af, ar = auroc(full), auroc(rest)

        rel = os.path.relpath(f, a.runs)
        if a.diagonal_only and nm not in rel:
            continue
        print(f'{rel:52s} {nm:20s} {jac:6.2f} {len(rows):7d} '
              f'{af:7.4f} {label_condition(nm, af):>9s} {ar:11.4f} '
              f'{IDENTITY[nm][1]:9.4f}')

    print('\nrestricted = homogeneous-protein sites removed (same subset as the '
          'identity baseline).\nid-restr = identity baseline on that same subset.')


if __name__ == '__main__':
    main()

