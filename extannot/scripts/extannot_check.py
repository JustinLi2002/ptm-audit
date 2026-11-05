"""External-annotation check: Step 1 (counts), Step 2 (rescoring of the deposited predictions) and the enrichment
analysis of PREREG.md (ptm-audit/extannot_20261102). No model is retrained.

ZENODO_X holds the extracted archives of the data record (rebuilt/, replica/, pdisjoint/, pdisjoint_runs_v2/);
SITES_TSV is the output of build_sites.py.

Usage: python -I extannot_check.py ZENODO_X SITES_TSV OUT_DIR"""
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ZX, SITES, OUT = map(Path, sys.argv[1:4])
OUT.mkdir(parents=True, exist_ok=True)
RUNS = ZX / "pdisjoint_runs_v2"

TASKS = ["phosphorylation_st", "phosphorylation_y", "acetylation_k", "methylation_k", "methylation_r",
         "sumoylation_k", "ubiquitination_k", "glycosylation_n"]
UNIT = {t: "methylation" if t.startswith("methylation") else t for t in TASKS}
KINDS = {"phosphorylation_st": {"phospho"}, "phosphorylation_y": {"phospho"}, "acetylation_k": {"acetyl"},
         "methylation_k": {"methyl_K", "methyl_R"}, "methylation_r": {"methyl_R"}, "sumoylation_k": {"sumo"},
         "ubiquitination_k": {"ubiquitin"}, "glycosylation_n": {"nglyc"}}
RESIDUES = {"phosphorylation_st": "ST", "phosphorylation_y": "Y", "acetylation_k": "K", "methylation_k": "KR",
            "methylation_r": "R", "sumoylation_k": "K", "ubiquitination_k": "K", "glycosylation_n": "N"}
# rebuild_datasets.py: the >=10-site donor threshold is evaluated on the merged modification type
THRESHOLD_GROUP = {"acetylation_k": ["acetylation_k"], "glycosylation_n": ["glycosylation_n"],
                   "methylation_k": ["methylation_k"], "methylation_r": ["methylation_k"],
                   "phosphorylation_st": ["phosphorylation_st", "phosphorylation_y"],
                   "phosphorylation_y": ["phosphorylation_st", "phosphorylation_y"],
                   "sumoylation_k": ["sumoylation_k"], "ubiquitination_k": ["ubiquitination_k"]}
MIN_SITES = 10
CONS = {"T": "replica", "U": "rebuilt"}
CELLS = [("TU", "replica", "rebuilt"), ("TT", "replica", "replica"), ("UU", "rebuilt", "rebuilt")]
FEATS = {"ESM-2": "__esm", "interaction": "", "ProtT5": "__prott5"}
SPLITS = (0, 1, 2)


def read_fasta(path):
    d, acc = {}, None
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if line.startswith(">"):
            acc = line[1:].split()[0]
            d[acc] = []
        elif acc:
            d[acc].append(line)
    return {k: "".join(v) for k, v in d.items()}


def auroc(y, s):
    y = np.asarray(y, bool)
    n1, n0 = y.sum(), (~y).sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(s)
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


ref = read_fasta(ZX / "pdisjoint" / "proteins.fasta")
report = [f"reference: {len(ref)} proteins (pdisjoint/proteins.fasta)"]

# ------------------------------------------------------------------ external sites, matched to our reference
sites = pd.read_csv(SITES, sep="\t", dtype={"window15": str}, keep_default_na=False)
match = Counter()
valid = defaultdict(set)  # (source set name) -> {(task, acc, pos)}
cplm_win = Counter()
for r in sites.itertuples(index=False):
    seq = ref.get(r.accession)
    if seq is None:
        match[(r.source, r.tier, "accession not in reference")] += 1
        continue
    if not (0 < r.position <= len(seq)) or seq[r.position - 1] != r.residue:
        match[(r.source, r.tier, "residue mismatch (dropped)")] += 1
        continue
    match[(r.source, r.tier, "matched")] += 1
    if r.source == "CPLM" and r.window15:
        p = r.position
        cplm_win["window agrees" if seq[p - 8:p + 7] == r.window15 else "window differs"] += 1
    sets = [f"{r.source}:{r.tier}"]
    if r.tier == "primary":
        sets.append("union")
    if r.source == "UniProt":
        sets.append("UniProt:primary+similarity")
    for t in TASKS:
        if r.kind in KINDS[t] and r.residue in RESIDUES[t]:
            for s in sets:
                valid[s].add((t, r.accession, r.position))
SOURCE_SETS = ["union"] + sorted(s for s in valid if s != "union")
report.append("external sites against the reference (residue check):")
report += [f"  {s}\t{t}\t{k}\t{n}" for (s, t, k), n in sorted(match.items())]
report.append(f"  CPLM, residue-matched sites, 15-residue window vs reference (diagnostic, post hoc): {dict(cplm_win)}")
report.append("source sets: " + ", ".join(SOURCE_SETS))
HIT = {s: {t: {(a, p) for (tt, a, p) in valid[s] if tt == t} for t in TASKS} for s in SOURCE_SETS}


def hits(df, task, source):
    h = HIT[source][task]
    return np.fromiter(((a, p) in h for a, p in zip(df.protein, df.pos)), bool, len(df))


# ------------------------------------------------------------------ Step 1: counts
data = {c: {t: pd.read_csv(ZX / c / f"{t}_all.tsv", sep="\t", usecols=["protein", "aa", "pos", "y"])
            for t in TASKS} for c in CONS.values()}
for t in TASKS:
    for c in CONS.values():
        d = data[c][t]
        bad = sum(ref[a][p - 1] != aa for a, aa, p in zip(d.protein, d.aa, d.pos))
        assert bad == 0, (c, t, bad)
    pr = data["replica"][t].query("y == 1")[["protein", "pos"]].sort_values(["protein", "pos"]).values
    pu = data["rebuilt"][t].query("y == 1")[["protein", "pos"]].sort_values(["protein", "pos"]).values
    assert (pr == pu).all(), t
report.append("checked: every row's residue agrees with the reference; positives identical in T and U")

donor = {}
for t in TASKS:
    cnt = Counter()
    for g in THRESHOLD_GROUP[t]:
        cnt.update(data["rebuilt"][g].query("y == 1").protein)
    donor[t] = {p for p, n in cnt.items() if n >= MIN_SITES}
    off = (~data["replica"][t].query("y == 0").protein.isin(donor[t])).sum()
    assert off == 0, (t, off)
report.append("checked: every threshold-sampled negative lies on a donor protein (>=10 merged-type positives)")

split_of = {s: pd.read_csv(ZX / "pdisjoint" / f"split_seed{s}.csv").set_index("accession")["split"] for s in SPLITS}
step1, donor_rows = [], []
for t in TASKS:
    for cl, c in CONS.items():
        d = data[c][t]
        flags = {s: hits(d, t, s) for s in SOURCE_SETS}
        for sp in SPLITS:
            part = d.protein.map(split_of[sp]).values
            assert not pd.isna(part).any()
            for pt in ("train", "val", "test"):
                m = part == pt
                neg, pos = m & (d.y.values == 0), m & (d.y.values == 1)
                for s in SOURCE_SETS:
                    step1.append(dict(task=t, split=sp, construction=cl, partition=pt, source=s,
                                      n_neg=int(neg.sum()), neg_hit=int((neg & flags[s]).sum()),
                                      n_pos=int(pos.sum()), pos_hit=int((pos & flags[s]).sum())))
                if cl == "U" and pt == "test":
                    isd = d.protein.isin(donor[t]).values
                    for s in SOURCE_SETS:
                        for lab, dm in (("donor", isd), ("non-donor", ~isd)):
                            donor_rows.append(dict(task=t, split=sp, source=s, protein_status=lab,
                                                   n_neg=int((neg & dm).sum()),
                                                   neg_hit=int((neg & dm & flags[s]).sum())))
step1 = pd.DataFrame(step1)
step1["neg_share"] = step1.neg_hit / step1.n_neg
step1["pos_share"] = step1.pos_hit / step1.n_pos
step1.to_csv(OUT / "step1_counts.tsv", sep="\t", index=False)
donor_df = pd.DataFrame(donor_rows)
donor_df["neg_share"] = donor_df.neg_hit / donor_df.n_neg
donor_df.to_csv(OUT / "step1_test_negatives_by_donor_status.tsv", sep="\t", index=False)


# ------------------------------------------------------------------ Step 2: rescoring, and enrichment
def pred(task, train, cond, sp, test, suffix=""):
    f = RUNS / f"{task}__{train}__{cond}__split{sp}{suffix}__on_{test}.pred.tsv"
    return pd.read_csv(f, sep="\t", usecols=["protein", "aa", "pos", "y", "y_pred"])


def run_json(task, train, cond, sp, suffix=""):
    return json.loads((RUNS / f"{task}__{train}__{cond}__split{sp}{suffix}.json").read_text(encoding="utf-8"))


long, enrich = [], []
max_dev = 0.0
for cell, train, test in CELLS:
    for t in TASKS:
        for sp in SPLITS:
            b = pred(t, train, "baseline", sp, test)
            key = list(zip(b.protein, b.pos))
            part = b.protein.map(split_of[sp])
            assert (part == "test").all(), (cell, t, sp)
            dt = data[test][t]
            dt = dt[dt.protein.map(split_of[sp]) == "test"]
            assert sorted(zip(dt.protein, dt.pos, dt.y)) == sorted(zip(b.protein, b.pos, b.y)), (cell, t, sp)
            scores = {"baseline": b.y_pred.values}
            for fn, suf in FEATS.items():
                f = pred(t, train, "ppi", sp, test, suf)
                assert list(zip(f.protein, f.pos)) == key and (f.y.values == b.y.values).all()
                scores[fn] = f.y_pred.values
            y = b.y.values == 1
            # the deposited AUROC (ensemble) must be reproduced before anything is corrected
            for fn, suf in [("baseline", "")] + list(FEATS.items()):
                js = run_json(t, train, "baseline" if fn == "baseline" else "ppi", sp, suf)
                max_dev = max(max_dev, abs(js["tests"][test]["auroc"] - auroc(y, scores[fn])))
            for s in [None] + SOURCE_SETS:
                h = hits(b, t, s) & ~y if s else np.zeros(len(b), bool)
                corrections = [("original", y, np.ones(len(b), bool))] if s is None else \
                    [("relabel", y | h, np.ones(len(b), bool)), ("drop", y, ~h)]
                for corr, yy, keep in corrections:
                    a0 = auroc(yy[keep], scores["baseline"][keep])
                    for fn in FEATS:
                        a1 = auroc(yy[keep], scores[fn][keep])
                        long.append(dict(cell=cell, task=t, split=sp, feature=fn, source=s or "none",
                                         correction=corr, n_test=int(keep.sum()), n_flagged=int(h.sum()),
                                         auroc_baseline=a0, auroc_feature=a1, delta=a1 - a0))
                if s and cell == "TU":
                    neg = ~y
                    row = dict(task=t, split=sp, source=s, n_neg=int(neg.sum()), n_flagged=int(h.sum()))
                    for fn in ["baseline"] + list(FEATS):
                        row[f"enrich_{fn}"] = auroc(h[neg], scores[fn][neg])
                    enrich.append(row)
report.append(f"checked: deposited ensemble AUROCs reproduced from y_pred, max |difference| = {max_dev:.2e}")
assert max_dev < 1e-4
long = pd.DataFrame(long)
long.to_csv(OUT / "step2_rescoring_by_split.tsv", sep="\t", index=False)
enrich = pd.DataFrame(enrich)
enrich.to_csv(OUT / "enrichment_by_split.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ summaries and the decision rule
per_task = long.groupby(["cell", "feature", "source", "correction", "task"], as_index=False).agg(
    delta=("delta", "mean"), n_flagged=("n_flagged", "sum"), n_test=("n_test", "sum"))
orig = per_task[per_task.correction == "original"].set_index(["cell", "feature", "task"]).delta
summ = []
for (cell, fn, s, corr), g in per_task.groupby(["cell", "feature", "source", "correction"]):
    g = g.set_index("task").loc[TASKS]
    unit_mean = g.delta.groupby(g.index.map(UNIT)).mean().mean()
    o = orig.loc[cell, fn].loc[TASKS]
    o_unit = o.groupby(o.index.map(UNIT)).mean().mean()
    summ.append(dict(cell=cell, feature=fn, source=s, correction=corr, mean8=g.delta.mean(), mean7=unit_mean,
                     tasks_negative=int((g.delta < 0).sum()), retained=unit_mean / o_unit,
                     n_flagged=int(g.n_flagged.sum()), **{t: g.delta[t] for t in TASKS}))
summ = pd.DataFrame(summ)
summ.to_csv(OUT / "step2_summary.tsv", sep="\t", index=False)

en = enrich.groupby(["source", "task"], as_index=False)[[c for c in enrich if c.startswith("enrich_")]].mean()
en["ESM-2_minus_baseline"] = en["enrich_ESM-2"] - en["enrich_baseline"]
en.to_csv(OUT / "enrichment_summary.tsv", sep="\t", index=False)

pri = summ.query("cell == 'TU' and feature == 'ESM-2' and source == 'union' and correction == 'relabel'").iloc[0]
if pri.mean8 >= 0:
    verdict = "EXPLAINED"
elif pri.retained >= 0.5 and pri.tasks_negative >= 6:
    verdict = "HARM SURVIVES"
else:
    verdict = "PARTLY EXPLAINED"
report.append("")
report.append("PRIMARY (TU, ESM-2, union of sources, relabel): "
              f"8-task mean {pri.mean8:+.4f}, 7-unit mean {pri.mean7:+.4f}, retained {pri.retained:.3f}, "
              f"{pri.tasks_negative}/8 tasks negative -> {verdict}")
(OUT / "report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
print("\n".join(report))
