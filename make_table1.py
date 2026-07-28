#!/usr/bin/env python3
"""
make_table1.py — assemble the main result table for paper 1.

Reads inference_out/summary.json (produced by rerun_inference.py) and combines it
with the identity-baseline numbers from the audit to produce:

  A. per-condition ensemble AUROC with across-seed SD
  B. the central comparison: how much the sequence model vs the identity baseline
     degrade when label-homogeneous proteins are removed
  C. shuffle control: fraction of the real feature's gain recovered by a
     randomly permuted vector
  D. paired test across seeds (baseline vs +PPI)

Usage:
    python make_table1.py --summary inference_out/summary.json
"""

import argparse
import json
import math
from collections import defaultdict

# identity baseline: (full test set, homogeneous proteins removed)
IDENTITY = {
    'acetylation_k': (0.8765, 0.7494), 'glycosylation_n': (0.9247, 0.5914),
    'methylation_k': (0.9571, 0.7830), 'methylation_r': (0.9450, 0.8001),
    'phosphorylation_st': (0.7729, 0.7300), 'phosphorylation_y': (0.8020, 0.7527),
    'sumoylation_k': (0.8201, 0.7092), 'ubiquitination_k': (0.7990, 0.7229),
}
# share of training sites on pure-positive proteins (audit table T3)
PUREPOS = {
    'acetylation_k': 22.5, 'glycosylation_n': 38.3, 'methylation_k': 25.2,
    'methylation_r': 29.8, 'phosphorylation_st': 4.4, 'phosphorylation_y': 2.4,
    'sumoylation_k': 17.9, 'ubiquitination_k': 13.5,
}
ORDER = ['glycosylation_n', 'methylation_r', 'methylation_k', 'acetylation_k',
         'sumoylation_k', 'ubiquitination_k', 'phosphorylation_st',
         'phosphorylation_y']
SHORT = {'glycosylation_n': 'N-Glyc N', 'methylation_r': 'Meth R',
         'methylation_k': 'Meth K', 'acetylation_k': 'Acetyl K',
         'sumoylation_k': 'Sumo K', 'ubiquitination_k': 'Ubiq K',
         'phosphorylation_st': 'Phospho S/T', 'phosphorylation_y': 'Phospho Y'}


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return num / den if den else float('nan')


def paired_t(xs, ys):
    """Two-sided paired t-test; returns (t, df). p looked up separately."""
    d = [x - y for x, y in zip(xs, ys)]
    n = len(d)
    if n < 2:
        return float('nan'), 0
    m = sum(d) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in d) / (n - 1))
    return (m / (sd / math.sqrt(n)) if sd else float('inf')), n - 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--summary', default='inference_out/summary.json')
    a = ap.parse_args()

    with open(a.summary) as fh:
        rows = json.load(fh)
    R = defaultdict(dict)
    for r in rows:
        R[r['ptm']][r['cond']] = r

    # ── A. per-condition AUROC ───────────────────────────────────────────────
    print('=== A. ensemble AUROC (across-seed SD in parentheses) ===')
    print(f"{'PTM':13s} {'baseline':>18s} {'+kinase':>18s} {'+PPI':>18s} "
          f"{'shuffled kin':>18s}")
    for p in ORDER:
        line = f'{SHORT[p]:13s}'
        for c in ('baseline', 'kinase', 'ppi', 'shuffled'):
            r = R[p].get(c)
            line += (f"  {r['ensemble']:.4f}({r['seed_sd']:.4f})" if r
                     else f"{'-':>18s}")
        print(line)

    # ── B. the central comparison ────────────────────────────────────────────
    print('\n=== B. degradation when label-homogeneous proteins are removed ===')
    print(f"{'PTM':13s} {'seq full':>9s} {'seq rest':>9s} {'seq drop':>9s} | "
          f"{'id full':>8s} {'id rest':>8s} {'id drop':>8s} | {'ratio':>7s} {'n':>7s}")
    sd_, idd_ = [], []
    for p in ORDER:
        r = R[p].get('baseline')
        if not r:
            continue
        sf, sr = r['ensemble'], r['restricted']
        idf, idr = IDENTITY[p]
        ds, di = sf - sr, idf - idr
        sd_.append(ds)
        idd_.append(di)
        ratio = di / ds if ds > 1e-9 else float('inf')
        print(f'{SHORT[p]:13s} {sf:9.4f} {sr:9.4f} {ds:+9.4f} | '
              f'{idf:8.4f} {idr:8.4f} {di:+8.4f} | {ratio:7.1f}x '
              f'{r["n_restricted"]:7d}')
    ms, mi = sum(sd_) / len(sd_), sum(idd_) / len(idd_)
    print(f'\n  mean drop: sequence {ms:+.4f}, identity {mi:+.4f}  '
          f'({mi / ms:.1f}x)')
    print('  the sequence model is stable on the restricted subset, so the subset\n'
          '  is not intrinsically harder; the identity baseline\'s drop is the artefact.')

    # dose-response
    xs = [PUREPOS[p] for p in ORDER if R[p].get('baseline')]
    ys = [IDENTITY[p][0] - IDENTITY[p][1] for p in ORDER if R[p].get('baseline')]
    print(f'\n  pure-positive share vs identity drop: r = {pearson(xs, ys):.3f} '
          f'(n={len(xs)})')

    # ── C. shuffle control ───────────────────────────────────────────────────
    print('\n=== C. shuffle control: gain recovered by a random vector ===')
    print(f"{'PTM':13s} {'baseline':>9s} {'real kin':>9s} {'shuffled':>9s} "
          f"{'recovered':>10s}")
    rec = []
    for p in ORDER:
        b, k, s = R[p].get('baseline'), R[p].get('kinase'), R[p].get('shuffled')
        if not (b and k and s):
            continue
        dr, ds = k['ensemble'] - b['ensemble'], s['ensemble'] - b['ensemble']
        pct = 100 * ds / dr if abs(dr) > 1e-9 else float('nan')
        rec.append(pct)
        print(f"{SHORT[p]:13s} {b['ensemble']:9.4f} {k['ensemble']:9.4f} "
              f"{s['ensemble']:9.4f} {pct:9.1f}%")
    if rec:
        print(f'\n  recovered {min(rec):.1f}%-{max(rec):.1f}% '
              f'(mean {sum(rec) / len(rec):.1f}%)')

    # ── D. paired test across seeds ──────────────────────────────────────────
    print('\n=== D. paired t across the 10 seeds (baseline vs +PPI) ===')
    print(f"{'PTM':13s} {'base mean':>10s} {'ppi mean':>10s} {'diff':>8s} "
          f"{'t':>8s} {'df':>4s}")
    for p in ORDER:
        b, q = R[p].get('baseline'), R[p].get('ppi')
        if not (b and q):
            continue
        xb, xq = b['seed_aurocs'], q['seed_aurocs']
        t, df = paired_t(xq, xb)
        print(f'{SHORT[p]:13s} {sum(xb) / len(xb):10.4f} {sum(xq) / len(xq):10.4f} '
              f'{sum(xq) / len(xq) - sum(xb) / len(xb):+8.4f} {t:8.2f} {df:4d}')
    print('\n  CAVEAT: seeds differ only in initialisation, not in the data split.\n'
          '  This tests consistency across initialisations, NOT generalisation to a\n'
          '  new split. Report it as such; the protein-disjoint runs are the real test.')


if __name__ == '__main__':
    main()
