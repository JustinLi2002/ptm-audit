#!/usr/bin/env python3
"""
False-negative sensitivity of the natural-negative result.

The harm reported under natural negatives is measured against labels that are
themselves incomplete: a curated site set records modifications that were
assayed, and assay coverage rises with annotation depth.  The same is true of
the relationship the model is said to invert -- the positive rate per protein
rises with depth partly because deeply annotated proteins were measured more.
This script asks how much of the reported harm survives if that circularity is
taken seriously.

Parameterisation
----------------
Let f(d) be the OBSERVED positive rate among candidate sites, as a function of
merged-type annotation depth d, and let f_ref be its value in the deepest
stratum.  Write the TRUE rate as

    g(d) = f_ref * (f(d) / f_ref) ** alpha
    pi(d) = f(d) / g(d) = (f(d) / f_ref) ** (1 - alpha)

where pi(d) is the probability that a true site at depth d was annotated.

    alpha = 1  ->  g = f, pi == 1.  No false negatives.  This reproduces the
                   published analysis exactly and is the script's self-check.
    alpha = 0  ->  g == f_ref everywhere.  The entire observed depth-rate
                   relationship is study bias; shallow proteins are as heavily
                   modified as deep ones, and all of the difference is missing
                   annotation.  This is the adversarial extreme.

The posterior that a site observed as negative is truly positive:

    q(d) = g(1 - pi) / [ g(1 - pi) + (1 - g) ]

Observed positives are never flipped (curated MS sites; false positives assumed
negligible -- state this in the manuscript).  Negatives in the rebuilt sets are
a uniform, label-blind subsample of candidates (rebuild_datasets.py:140), so
q(d) applies to them unchanged; the denominator of f(d), however, must be the
full candidate count, not the retained negative count.

Usage
-----
    python fn_sensitivity.py --selfcheck        # verify wiring, run nothing
    python fn_sensitivity.py --arm esm
    python fn_sensitivity.py --arm ppi --draws 200 --out ../results/
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

# ---------------------------------------------------------------------------
# adapter.  These are the only project-specific bindings; everything below is
# generic.  Definitions are imported rather than restated so that d matches
# Figure 3 exactly.
# ---------------------------------------------------------------------------
BASE = Path("/home/FCAM/juli/HRP")
RUNS = BASE / "pdisjoint_runs_v2"
sys.path.insert(0, str(BASE / "ptm-audit" / "analysis"))
sys.path.insert(0, str(BASE / "ptm-audit"))

PTMS = MERGED = depths = rb = None


def load_project_modules(arm):
    """Import the project's own definitions rather than restating them.

    mechanism_chain parses sys.argv at module level (FEAT is a module global),
    so importing it with our own flags present makes ITS argparse fail.  Swap
    argv for the duration of the import.
    """
    global PTMS, MERGED, depths, rb
    saved = sys.argv
    sys.argv = [saved[0], "--feat", arm]
    try:
        from mechanism_chain import PTMS as _P, MERGED as _M, depths as _d
        import rebuild_datasets as _rb
    finally:
        sys.argv = saved
    PTMS, MERGED, depths, rb = _P, _M, _d, _rb

SPLITS = (0, 1, 2)
TRAIN = "replica"      # threshold-sampled training
EVAL = "rebuilt"       # natural negatives


def pred_path(task, cond, split, arm):
    """cond is 'baseline' or 'ppi'; arm selects the feature family."""
    suffix = "__esm" if (arm == "esm" and cond != "baseline") else ""
    return RUNS / f"{task}__{TRAIN}__{cond}__split{split}{suffix}__on_{EVAL}.pred.tsv"


_FASTA_CACHE = {}


def load_fasta(known_ids):
    """{protein: sequence}, with the header ID convention inferred, not guessed.

    rebuild_datasets exposes the path but builds the dict in its own __main__,
    so we parse it here.  UniProt FASTA headers appear both as bare accessions
    and as sp|ACC|NAME; rather than assume, all plausible extractions are tried
    and the one that best covers the proteins actually present in the
    predictions is kept.  Coverage is reported so a bad guess cannot pass
    silently.
    """
    if "fa" in _FASTA_CACHE:
        return _FASTA_CACHE["fa"]

    path = None
    for attr in ("FASTA", "FASTA_PATH", "PROTEOME", "SEQS"):
        v = getattr(rb, attr, None)
        if isinstance(v, (str, Path)) and Path(v).exists():
            path = Path(v)
            break
    if path is None:
        raise RuntimeError("could not locate the proteome FASTA on "
                           "rebuild_datasets; pass --fasta")

    headers, seqs = [], []
    with open(path) as fh:
        cur = []
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
                headers.append(line[1:].strip())
            else:
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur))
    if len(headers) != len(seqs):
        raise RuntimeError(f"{path}: {len(headers)} headers, {len(seqs)} seqs")

    def first_token(h):
        return h.split()[0] if h.split() else ""

    schemes = {
        "first token": first_token,
        "pipe field 1": lambda h: (first_token(h).split("|")[1]
                                   if first_token(h).count("|") >= 2 else ""),
        "pipe field 0": lambda h: first_token(h).split("|")[0],
        "before dash": lambda h: first_token(h).split("-")[0],
    }
    known = set(known_ids)
    best, best_cov, report = None, -1, []
    for name, fn in schemes.items():
        ids = [fn(h) for h in headers]
        cov = len(known & set(ids))
        report.append((name, cov))
        if cov > best_cov:
            best, best_cov, best_ids = name, cov, ids

    _FASTA_CACHE["report"] = (path, len(headers), report, best, best_cov,
                              len(known))
    if best_cov < 0.95 * len(known):
        raise RuntimeError(
            f"no header convention covers the predicted proteins "
            f"({best_cov}/{len(known)} with '{best}'); tried {report}")

    fa = {}
    for i, s in zip(best_ids, seqs):
        if i and i not in fa:
            fa[i] = s
    _FASTA_CACHE["fa"] = fa
    return fa


def candidate_counts(task, proteins):
    """Candidate sites per protein, from the same enumeration the rebuild used.

    Returns {protein: n_candidates}.  Includes the positives, since f(d) is a
    rate over candidates.
    """
    fa = load_fasta(proteins)
    cands = rb.enumerate_candidates(fa, task, set(proteins))
    n = defaultdict(int)
    for c in cands:
        n[c[0]] += 1
    return n


# ---------------------------------------------------------------------------
# core
# ---------------------------------------------------------------------------
def depth_bins(depth_of, proteins, n_bins):
    """Quantile bins over merged-type depth, weighted by protein count.

    Depth in the natural-negative partitions runs 1..~15 with a low median, so
    equal-width bins would leave the deep strata nearly empty.  Bins are merged
    where quantiles collide on ties.
    """
    d = np.array([depth_of.get(p, 0) for p in proteins], float)
    edges = np.unique(np.quantile(d, np.linspace(0, 1, n_bins + 1)))
    if len(edges) < 3:
        raise ValueError("depth has too few distinct values to bin")
    return edges, np.clip(np.digitize(d, edges[1:-1], right=True),
                          0, len(edges) - 2)


def observed_rates(frame, cand_of, depth_of, n_bins):
    """f(b) per depth bin: annotated positives / candidate sites.

    The candidate denominator matters: the retained negatives are a subsample,
    so positives/(positives+retained negatives) would overstate f and shrink
    the estimated false-negative burden -- i.e. bias the analysis in the
    direction we are trying to defend against.
    """
    per_prot = frame.groupby("protein")["y"].sum()
    proteins = per_prot.index.to_numpy()
    edges, b = depth_bins(depth_of, proteins, n_bins)

    n_bin = len(edges) - 1
    pos = np.zeros(n_bin)
    cand = np.zeros(n_bin)
    for p, bi in zip(proteins, b):
        pos[bi] += per_prot[p]
        cand[bi] += cand_of.get(p, 0)
    if (cand == 0).any():
        raise ValueError("a depth bin has no candidate sites; check the adapter")

    f = pos / cand
    bin_of = dict(zip(proteins, b))
    return f, bin_of, edges


def posterior(f, alpha):
    """q per bin.

    The reference stratum is the one with the HIGHEST observed rate, not the
    deepest.  Observed rate per candidate site is not monotone in merged
    annotation depth -- it peaks mid-range and falls in the deepest strata,
    because candidate count grows with protein length faster than annotated
    sites do in heavily studied proteins.  Anchoring on the deepest bin would
    therefore require clipping ratios above one, which silently sets q to zero
    in exactly the strata that clip.

    Anchoring on the maximum needs no monotonicity, needs no clipping, and is
    the most adversarial choice available in this family: every other anchor
    implies a smaller false-negative burden.
    """
    ref = int(np.argmax(f))
    f_ref = f[ref]
    ratio = f / f_ref
    if ratio.max() > 1 + 1e-12:
        raise AssertionError("ratio exceeds one after anchoring on the max")
    g = f_ref * ratio ** alpha
    pi = ratio ** (1.0 - alpha)
    q = g * (1 - pi) / (g * (1 - pi) + (1 - g))
    return q, ref


def delta_under_flips(frame, q_of_site, rng, draws):
    """Monte Carlo Delta AUROC / Delta AUPRC over resampled labels."""
    y = frame["y"].to_numpy()
    base = frame["score_base"].to_numpy()
    aug = frame["score_aug"].to_numpy()
    neg = y == 0

    d_auroc = np.empty(draws)
    d_auprc = np.empty(draws)
    for i in range(draws):
        yy = y.copy()
        flip = rng.random(neg.sum()) < q_of_site[neg]
        yy[np.flatnonzero(neg)[flip]] = 1
        if yy.all() or not yy.any():
            d_auroc[i] = d_auprc[i] = np.nan
            continue
        d_auroc[i] = roc_auc_score(yy, aug) - roc_auc_score(yy, base)
        d_auprc[i] = (average_precision_score(yy, aug)
                      - average_precision_score(yy, base))
    return d_auroc, d_auprc


def solve_alpha(f, bins, y, target_flipped):
    """alpha such that the expected number of imputed positives hits a target.

    A fixed alpha is not a fixed strength of doubt.  alpha acts on the dynamic
    range of f, which spans 1.5x in phosphorylation Y and 11x in
    phosphorylation S/T, so the same alpha implies a 1.14x inflation of the
    true site count in one task and 2.02x in another.  Fixing the inflation
    instead makes the assumption comparable across tasks and legible to a
    reader who has an opinion about how complete a given PTM's annotation is.

    exp_flipped is decreasing in alpha (q -> 0 as alpha -> 1), so bisect.
    """
    neg = y == 0

    def flipped(a):
        q, _ = posterior(f, a)
        return q[bins][neg].sum()

    hi_burden = flipped(0.0)
    if target_flipped >= hi_burden:
        return 0.0, float(hi_burden), False        # unreachable; report actual
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if flipped(mid) > target_flipped:
            lo = mid
        else:
            hi = mid
    a = (lo + hi) / 2
    return a, float(flipped(a)), True


def load_split(task, split, arm):
    b = pd.read_csv(pred_path(task, "baseline", split, arm), sep="\t",
                    usecols=["protein", "pos", "y", "y_pred"])
    a = pd.read_csv(pred_path(task, "ppi", split, arm), sep="\t",
                    usecols=["protein", "pos", "y", "y_pred"])
    m = b.merge(a, on=["protein", "pos"], suffixes=("_base", "_aug"))
    if len(m) != len(b) or len(m) != len(a):
        raise ValueError(f"{task} split{split}: site sets differ between "
                         f"baseline ({len(b)}) and augmented ({len(a)}); "
                         f"merged to {len(m)}")
    if not (m["y_base"] == m["y_aug"]).all():
        raise ValueError(f"{task} split{split}: labels disagree between files")
    return m.rename(columns={"y_base": "y",
                             "y_pred_base": "score_base",
                             "y_pred_aug": "score_aug"})


def run(arm, alphas, draws, n_bins, seed, out, inflations=None):
    _, merged_depth = depths()
    rows = []

    for task in PTMS:
        depth_of = merged_depth[MERGED[task]]
        for split in SPLITS:
            frame = load_split(task, split, arm)
            cand_of = candidate_counts(task, frame["protein"].unique())
            f, bin_of, edges = observed_rates(frame, cand_of, depth_of, n_bins)
            bins = frame["protein"].map(bin_of).to_numpy()

            grid = [(a, None) for a in alphas]
            n_pos = int((frame["y"] == 1).sum())
            _, ceiling, _ = solve_alpha(f, bins, frame["y"].to_numpy(), np.inf)
            max_infl = round(1.0 + ceiling / max(n_pos, 1), 4)
            if inflations:
                for lam in inflations:
                    a, got, ok = solve_alpha(f, bins, frame["y"].to_numpy(),
                                             (lam - 1) * n_pos)
                    grid.append((a, (lam, got / max(n_pos, 1) + 1.0, ok)))

            for alpha, infl in grid:
                q, ref_bin = posterior(f, alpha)
                rng = np.random.default_rng(
                    abs(hash((task, split, round(alpha, 4), seed))) % 2**32)
                da, dp = delta_under_flips(frame, q[bins], rng, draws)
                rows.append(dict(
                    task=task, split=split, alpha=alpha,
                    mode="alpha" if infl is None else "inflation",
                    target_inflation=np.nan if infl is None else infl[0],
                    actual_inflation=(np.nan if infl is None
                                      else round(infl[1], 4)),
                    inflation_reached=(True if infl is None else infl[2]),
                    max_inflation=max_infl,
                    n_sites=len(frame), n_proteins=frame["protein"].nunique(),
                    exp_flipped=float(q[bins][frame["y"] == 0].sum()),
                    d_auroc=np.nanmean(da), d_auprc=np.nanmean(dp),
                    d_auroc_lo=np.nanpercentile(da, 2.5),
                    d_auroc_hi=np.nanpercentile(da, 97.5),
                    d_auprc_lo=np.nanpercentile(dp, 2.5),
                    d_auprc_hi=np.nanpercentile(dp, 97.5),
                    f_shallow=float(f[0]), f_deep=float(f[-1]),
                    f_ref=float(f[ref_bin]), ref_bin=ref_bin,
                    f_monotone=bool(np.all(np.diff(f) >= 0)),
                ))
            print(f"  {task:20s} split{split}  "
                  f"n={len(frame):7d}  bins={len(f)}  "
                  f"f {f[0]:.4f}..{f[-1]:.4f}", flush=True)

    df = pd.DataFrame(rows)
    grid_rows = df[df["mode"] == "alpha"]
    by_alpha = (grid_rows.groupby("alpha")[["d_auroc", "d_auprc"]]
                  .mean().reset_index())

    print("\n=== mean over 8 tasks x 3 splits (alpha grid) ===")
    print(f"{'alpha':>6s} {'dAUROC':>9s} {'dAUPRC':>9s} {'tasks<0':>8s}")
    for a in alphas:
        s = grid_rows[grid_rows.alpha == a].groupby("task")["d_auroc"].mean()
        r = by_alpha[by_alpha.alpha == a].iloc[0]
        print(f"{a:6.2f} {r.d_auroc:+9.4f} {r.d_auprc:+9.4f} "
              f"{int((s < 0).sum()):5d}/8")

    if inflations:
        print("\n=== by target inflation (alpha solved per task/split) ===")
        print(f"{'lambda':>7s} {'alpha':>7s} {'actual':>7s} "
              f"{'dAUROC':>9s} {'dAUPRC':>9s} {'tasks<0':>8s}")
        ceil = (df.groupby("task")["max_inflation"].mean()
                  .sort_values())
        print("  per-task ceiling at alpha=0: "
              + ", ".join(f"{t.split('_')[0][:6]} {v:.2f}x"
                          for t, v in ceil.items()))
        print(f"  universally reachable target: <= {ceil.min():.2f}x\n")
        for lam in inflations:
            sub = df[(df["mode"] == "inflation")
                     & (df.target_inflation == lam)]
            s7 = sub.groupby("task")["d_auroc"].mean()
            print(f"{lam:7.2f} {sub.alpha.mean():7.3f} "
                  f"{sub.actual_inflation.mean():7.3f} "
                  f"{sub.d_auroc.mean():+9.4f} {sub.d_auprc.mean():+9.4f} "
                  f"{int((s7 < 0).sum()):5d}/8")
            capped = sorted(set(sub[sub.inflation_reached == False].task))
            if capped:
                print(f"        capped at alpha=0 (ceiling below target): "
                      f"{', '.join(c.split('_')[0][:6] for c in capped)}")

    if out:
        out = Path(out)
        out.mkdir(parents=True, exist_ok=True)
        df.to_csv(out / f"fn_sensitivity_{arm}.tsv", sep="\t", index=False)
        by_alpha.to_csv(out / f"fn_sensitivity_{arm}_summary.tsv",
                        sep="\t", index=False)
        print(f"\nwrote {out}/fn_sensitivity_{arm}.tsv")
    return df


def selfcheck(arm):
    """Verify the adapter before committing to a full run."""
    print("PTMS      :", list(PTMS))
    print("MERGED    :", dict(MERGED))
    g, m = depths()
    print(f"depths()  : {len(g)} proteins global, "
          f"{len(m)} merged groups {list(m)}")

    task = PTMS[0]
    p = pred_path(task, "baseline", 0, arm)
    print(f"pred path : {p}\n            exists={p.exists()}")
    p2 = pred_path(task, "ppi", 0, arm)
    print(f"augmented : {p2}\n            exists={p2.exists()}")
    if not (p.exists() and p2.exists()):
        print("!! path template wrong -- fix pred_path() before running")
        return

    frame = load_split(task, 0, arm)
    print(f"merged    : {len(frame)} sites, "
          f"{frame['protein'].nunique()} proteins, "
          f"pos rate {frame['y'].mean():.4f}")

    cand_of = candidate_counts(task, frame["protein"].unique())
    path, nh, report, best, cov, nk = _FASTA_CACHE["report"]
    print(f"fasta     : {path}  ({nh} records)")
    print(f"            coverage by scheme {report}")
    print(f"            using '{best}': {cov}/{nk} predicted proteins matched")
    obs = frame.groupby("protein").size()
    print(f"candidates: median {np.median(list(cand_of.values())):.0f} per "
          f"protein vs {obs.median():.0f} rows retained "
          f"(candidates must be >= rows)")
    bad = [p for p in obs.index if cand_of.get(p, 0) < obs[p]]
    print(f"            {len(bad)} proteins violate that -- must be 0")

    depth_of = m[MERGED[task]]
    d = np.array([depth_of.get(p, 0) for p in frame["protein"].unique()])
    print(f"depth     : min {d.min()} median {np.median(d):.0f} max {d.max()}, "
          f"{(d == 0).sum()} proteins at zero")

    f, bin_of, edges = observed_rates(frame, cand_of, depth_of, 6)
    b = np.array([bin_of[p] for p in frame["protein"].unique()])
    cand = np.array([cand_of.get(p, 0) for p in frame["protein"].unique()])
    print(f"f(d)      : {np.round(f, 4)} over edges {np.round(edges, 1)}")
    print(f"            median candidates per bin "
          f"{[int(np.median(cand[b == i])) for i in range(len(f))]}")
    print(f"            monotone increasing: {bool(np.all(np.diff(f) >= 0))} "
          f"(not required; anchor is the max)")
    q0, ref = posterior(f, 0.0)
    q1, _ = posterior(f, 1.0)
    print(f"reference : bin {ref} (f = {f[ref]:.4f})"
          + ("   !! max is the shallowest bin -- the depth-rate relationship "
             "runs backwards here, read the alpha=0 extreme with care"
             if ref == 0 else ""))
    print(f"q at a=0  : {np.round(q0, 4)}")
    print(f"q at a=1  : {np.round(q1, 4)}  (must be all zero)")
    nneg = (frame["y"] == 0).sum()
    exp = float(q0[frame["protein"].map(bin_of).to_numpy()][frame["y"] == 0].sum())
    print(f"burden    : at alpha=0, {exp:.0f} of {nneg} negatives flip "
          f"({exp / nneg:.1%}); observed positives {(frame['y'] == 1).sum()}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", choices=["esm", "ppi"], default="esm")
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--bins", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--alphas", type=str, default="0,0.25,0.5,0.75,1")
    ap.add_argument("--inflations", type=str, default=None,
                    help="comma-separated true/observed site-count ratios, "
                         "e.g. 1.25,1.5,2.0; alpha is solved per task")
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--selfcheck", action="store_true")
    a = ap.parse_args()
    load_project_modules(a.arm)

    if a.selfcheck:
        selfcheck(a.arm)
    else:
        run(a.arm, [float(x) for x in a.alphas.split(",")],
            a.draws, a.bins, a.seed, a.out,
            [float(x) for x in a.inflations.split(",")] if a.inflations
            else None)
