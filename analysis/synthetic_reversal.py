#!/usr/bin/env python3
"""Synthetic demonstration that the sign reversal is a property of the
construction, not of this biology.

The mechanism reported in the manuscript requires only three things, none of
which is specific to proteins or to post-translational modification:

  (i)   labels are defined on items that are grouped into entities;
  (ii)  negatives are constructed rather than observed, by a rule that draws
        them only from entities carrying at least T positives;
  (iii) a feature that is constant within an entity partially encodes how many
        positives that entity has.

This script generates data satisfying those three conditions and nothing else,
splits it disjointly by entity, and reproduces the full pattern: the
entity-constant channel appears beneficial when evaluated against
threshold-drawn negatives, harms performance against naturally drawn ones,
and is there outperformed by a permuted vector of the same shape.

Because T and the strength of (iii) are free parameters here, the script also
maps where in that plane the sign of the contribution changes, which the real
data cannot show: the real benchmark offers a single point in it.

Nothing in the generative model refers to sequence, structure or biology. The
generative parameters are calibrated to the values measured in the manuscript
(per-entity positive rate correlating with depth at rho ~ +0.7; feature-to-log-
depth ridge R^2 in 0.13-0.58; class ratio 1:5 without the threshold), but the
qualitative result does not depend on that calibration, which is the point of
the sweep.

    python synthetic_reversal.py --mode demo    # one configuration, ~1 min
    python synthetic_reversal.py --mode sweep   # T x R^2 grid, ~40 min on CPU

Outputs a TSV per mode and, for the sweep, a heatmap.
"""
import argparse
import itertools

import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

# ─────────────────────────────── generation ───────────────────────────────

def generate(n_entities=6000, feat_dim=16, local_dim=16, r2_target=0.30,
             prop_scale=1.15, rng=None):
    """One synthetic population.

    Each entity carries a latent character z that determines both how prone its
    items are to being positive and, through that, how many positives it ends up
    with. The entity-constant feature encodes z with a controllable fidelity.
    This is the structure the real case has: an embedding reflects a protein's
    biology, which drives both its modification propensity and how many sites
    have been found on it. The feature is therefore genuinely informative --
    which is why it helps once the sampling is corrected -- while also carrying
    the quantity the threshold rule turns into a label.

    Annotation depth is the number of positives an entity happens to have: an
    observed consequence, not an input.
    """
    rng = rng or np.random.default_rng(0)

    n_items = 5 + rng.negative_binomial(3, 0.12, n_entities)
    z = rng.normal(0, 1, n_entities)                     # entity character
    prop = expit(-1.6 + prop_scale * z)                  # true propensity

    w = rng.normal(0, 1, local_dim) / np.sqrt(local_dim)

    ent, x_local, y = [], [], []
    for e in range(n_entities):
        m = n_items[e]
        xi = rng.normal(0, 1, (m, local_dim))
        p = expit(xi @ w * 1.2 + logit(prop[e]))
        ent.append(np.full(m, e))
        x_local.append(xi)
        y.append(rng.random(m) < p)
    ent = np.concatenate(ent)
    x_local = np.vstack(x_local)
    y = np.concatenate(y).astype(int)

    depth = np.bincount(ent[y == 1], minlength=n_entities)

    # entity-constant feature: one informative direction carrying z, plus
    # isotropic noise. r2_target sets the fidelity; the R^2 actually achieved
    # against log depth is measured afterwards rather than assumed.
    direction = rng.normal(0, 1, feat_dim)
    direction /= np.linalg.norm(direction)
    feat = (np.sqrt(r2_target) * np.outer(z, direction)
            + np.sqrt(1 - r2_target) * rng.normal(0, 1, (n_entities, feat_dim)))

    return dict(ent=ent, x_local=x_local, y=y, depth=depth, feat=feat,
                n_entities=n_entities)


def measured_r2(feat, depth, seed=0):
    """Ridge R^2 from the entity feature to log depth, as in the manuscript."""
    pipe = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 13)))
    t = np.log1p(depth)
    p = cross_val_predict(pipe, feat, t, cv=KFold(5, shuffle=True,
                                                  random_state=seed))
    return 1 - ((t - p) ** 2).sum() / ((t - t.mean()) ** 2).sum()


# ─────────────────────────────── construction ─────────────────────────────

def construct(d, threshold, ratio=5, rng=None):
    """Build a dataset by the negative-sampling rule under study.

    threshold = 0 draws negatives uniformly from every entity (the natural
    construction). threshold = T draws them only from entities carrying at
    least T positives, so entities below T contribute positives and no
    negatives, and every one of their items therefore carries the same label.
    """
    rng = rng or np.random.default_rng(1)
    pos = np.where(d['y'] == 1)[0]
    neg_pool = np.where(d['y'] == 0)[0]
    if threshold > 0:
        donor = d['depth'] >= threshold
        neg_pool = neg_pool[donor[d['ent'][neg_pool]]]
    k = min(len(neg_pool), ratio * len(pos))
    neg = rng.choice(neg_pool, k, replace=False)
    idx = np.concatenate([pos, neg])
    rng.shuffle(idx)
    return idx


def homogeneity(d, idx):
    """Share of selected items lying on entities that are pure-positive."""
    df = pd.DataFrame({'e': d['ent'][idx], 'y': d['y'][idx]})
    g = df.groupby('e')['y'].agg(['mean', 'size'])
    pure = g[g['mean'] == 1.0]['size'].sum()
    return 100 * pure / len(idx)


# ─────────────────────────────── evaluation ───────────────────────────────

def fit_predict(d, tr_idx, te_idx, mode, rng=None, seed=0):
    """mode: 'base' (local only), 'real', or 'perm'."""
    rng = rng or np.random.default_rng(seed)
    f = d['feat']
    if mode == 'perm':
        f = f[rng.permutation(len(f))]

    def X(idx):
        if mode == 'base':
            return d['x_local'][idx]
        return np.hstack([d['x_local'][idx], f[d['ent'][idx]]])

    clf = make_pipeline(StandardScaler(),
                        MLPClassifier(hidden_layer_sizes=(64, 32),
                                      alpha=1e-3, max_iter=400,
                                      early_stopping=True, n_iter_no_change=12,
                                      random_state=seed))
    clf.fit(X(tr_idx), d['y'][tr_idx])
    return clf.predict_proba(X(te_idx))[:, 1]


def one_configuration(threshold, r2_target, seed=0, n_entities=6000):
    rng = np.random.default_rng(seed)
    d = generate(n_entities=n_entities, r2_target=r2_target, rng=rng)

    # entity-disjoint split
    perm = rng.permutation(d['n_entities'])
    test_ent = np.zeros(d['n_entities'], bool)
    test_ent[perm[:d['n_entities'] // 5]] = True
    in_test = test_ent[d['ent']]

    out = {'threshold': threshold, 'r2_target': r2_target, 'seed': seed,
           'r2_measured': measured_r2(d['feat'], d['depth'], seed)}

    # training set under the rule being tested
    sel_thr = construct(d, threshold, rng=rng)
    tr = sel_thr[~in_test[sel_thr]]
    out['homog_train'] = homogeneity(d, tr)

    # two test partitions over the SAME held-out entities, differing only in
    # how their negatives were drawn
    tests = {}
    for name, T in (('threshold', threshold), ('natural', 0)):
        sel = construct(d, T, rng=np.random.default_rng(seed + 99))
        tests[name] = sel[in_test[sel]]
        out[f'homog_test_{name}'] = homogeneity(d, tests[name])

    scores = {}
    for mode in ('base', 'real', 'perm'):
        for name, te in tests.items():
            p = fit_predict(d, tr, te, mode, seed=seed)
            scores[(mode, name)] = p
            out[f'auroc_{mode}_{name}'] = roc_auc_score(d['y'][te], p)
            out[f'auprc_{mode}_{name}'] = average_precision_score(d['y'][te], p)

    for name in tests:
        b = out[f'auroc_base_{name}']
        out[f'd_real_{name}'] = out[f'auroc_real_{name}'] - b
        out[f'd_perm_{name}'] = out[f'auroc_perm_{name}'] - b
        out[f'real_beats_perm_{name}'] = int(
            out[f'auroc_real_{name}'] > out[f'auroc_perm_{name}'])

    # does the channel invert the score-versus-depth relationship?
    te = tests['natural']
    df = pd.DataFrame({'e': d['ent'][te]})
    for mode in ('base', 'real'):
        df['s'] = scores[(mode, 'natural')]
        g = df.groupby('e')['s'].mean()
        out[f'rho_score_depth_{mode}'] = spearmanr(
            d['depth'][g.index], g.values).statistic
    out['rho_delta'] = out['rho_score_depth_real'] - out['rho_score_depth_base']

    # the true relationship, for reference
    df2 = pd.DataFrame({'e': d['ent'], 'y': d['y']})
    g2 = df2.groupby('e')['y'].mean()
    out['rho_truth'] = spearmanr(d['depth'][g2.index], g2.values).statistic
    return out


# ─────────────────────────────── drivers ──────────────────────────────────

def demo(args):
    rows = [one_configuration(t, 0.30, seed=s, n_entities=args.n_entities)
            for t in (0, 10) for s in range(args.seeds)]
    df = pd.DataFrame(rows)
    df.to_csv(f'{args.out}/synthetic_demo.tsv', sep='\t', index=False)
    m = df.groupby('threshold').mean(numeric_only=True)
    print('\nmeasured feature->log depth R^2 : %.3f' % df.r2_measured.mean())
    print('true rho(depth, positive rate)  : %+.3f' % df.rho_truth.mean())
    print('\n%-28s %12s %12s' % ('', 'no threshold', 'threshold 10'))
    for k, lab in [('homog_train', 'pure-positive train items %'),
                   ('homog_test_natural', 'pure-positive natural test %'),
                   ('homog_test_threshold', 'pure-positive thr. test %'),
                   ('d_real_threshold', 'dAUROC, threshold test'),
                   ('d_real_natural', 'dAUROC, natural test'),
                   ('d_perm_natural', 'dAUROC permuted, natural'),
                   ('real_beats_perm_natural', 'real > perm on natural'),
                   ('rho_score_depth_base', 'rho(score, depth) base'),
                   ('rho_score_depth_real', 'rho(score, depth) +feature')]:
        print('%-28s %12.3f %12.3f' % (lab, m.loc[0, k], m.loc[10, k]))
    print('\nWritten to %s/synthetic_demo.tsv' % args.out)


def sweep(args):
    grid = list(itertools.product((0, 1, 2, 5, 10, 20),
                                  (0.02, 0.05, 0.10, 0.20, 0.40, 0.60)))
    rows = []
    for i, (t, r) in enumerate(grid, 1):
        for s in range(args.seeds):
            rows.append(one_configuration(t, r, seed=s,
                                          n_entities=args.n_entities))
        print('  %d/%d  threshold=%2d  r2=%.2f' % (i, len(grid), t, r),
              flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(f'{args.out}/figure6.tsv', sep='\t', index=False)

    piv = df.pivot_table(index='threshold', columns='r2_target',
                         values='d_real_natural')
    print('\ndAUROC on the natural test partition:\n')
    print(piv.round(3).to_string())

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.4))
        for k, (val, title) in enumerate([
                ('d_real_natural', 'a  $\\Delta$AUROC, natural negatives'),
                ('rho_delta', 'b  change in $\\rho$(score, depth)')]):
            p = df.pivot_table(index='threshold', columns='r2_target', values=val)
            v = np.abs(p.values).max()
            im = ax[k].imshow(p.values, cmap='RdBu_r', vmin=-v, vmax=v,
                              aspect='auto', origin='lower')
            ax[k].set_xticks(range(len(p.columns)))
            meas = df.groupby('r2_target')['r2_measured'].mean()
            ax[k].set_xticklabels([f'{meas[c]:.2f}' for c in p.columns],
                                  fontsize=7)
            ax[k].set_yticks(range(len(p.index)))
            ax[k].set_yticklabels(p.index, fontsize=7)
            ax[k].set_xlabel('feature $\\to$ log depth, measured $R^2$', fontsize=8)
            ax[k].set_ylabel('donor threshold', fontsize=8)
            ax[k].set_title(title, loc='left', fontsize=8)
            fig.colorbar(im, ax=ax[k], fraction=0.046)
        fig.tight_layout()
        fig.savefig(f'{args.out}/figure6.png', dpi=300,
                    bbox_inches='tight')
        print('\nHeatmap written to %s/figure6.png' % args.out)
    except ImportError:
        pass


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', default='demo', choices=['demo', 'sweep'])
    ap.add_argument('--seeds', type=int, default=3)
    ap.add_argument('--n-entities', type=int, default=6000)
    ap.add_argument('--out', default='.')
    a = ap.parse_args()
    (demo if a.mode == 'demo' else sweep)(a)
