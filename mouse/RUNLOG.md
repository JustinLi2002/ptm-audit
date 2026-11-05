# Mouse replication: run log

## 2026-11-02

**Approval.** Justin approved PREREG_mouse.md as written ("都同意，下载吧，HRP 也跑"). The header now reads FIXED. Its sha256 is in PREREG_mouse.sha256, taken before the header edit; the approved text is otherwise unchanged.

**Reference download.** UniProtKB/Swiss-Prot *Mus musculus*, canonical, release 2026_03 (02-Sep-2026), from rest.uniprot.org. 17,283 entries; sha256 values are in `data/uniprot/sha256.txt`. The first, uncompressed request returned "Internal server error" and was discarded; the compressed request succeeded.

**Mapping** (`scripts/map_sites.py` → `build/mapping_report.json`):
- 25,140 of 29,054 prefiltered S/T sites map (86.5%, above the 80% stop rule);
- 7,069 S/T proteins, of which 490 have merged S/T/Y depth ≥ 10;
- 17,103 S/T sites lie on proteins below the threshold;
- 335 Y sites.

**Constructions.** `rebuild_datasets.py` (repo commit 8393691, unchanged), run locally:
- U: 25,140 positives + 125,700 negatives; positive rate 0.167; pure-positive sites 0.1%.
- T: 90,863 negatives (pool-limited; maximum ratio 3.6); positive rate 0.217; pure-positive sites 14.7%.

These are recorded in `build/construction_prediction.txt` before any split. The negatives depend on Python's string-hash order (set iteration), so the uploaded files *are* the construction. Regenerating them needs a fixed PYTHONHASHSEED or the files themselves.

**HRP submission** (`build/hrp/submit_mouse_v1.sh`, tarball sha256 4ec88f30…, unpacked to `/home/FCAM/juli/HRP/mouse_v1`):
- First connection: no output at all, and evidently never reached the remote, because `mouse_v1` did not exist on the second try.
- Second connection: rc 0, 49 s. Jobs:
  - prep 4437956 (CD-HIT, splits 0–2, matched);
  - esm 4437957;
  - train 4437958 (array 0–26%9, afterok on both).
- Repository scripts on HRP match the local LF hashes:
  - train_pdisjoint.py f381a95b…
  - make_pdisjoint_split.py 6dfe95b3…
  - extract_esm.py 31586bbe…
- ESM-2 650M weights are in the torch hub cache.

**Analysis script** (`scripts/analyse_mouse.py`) was written before any result. A dry run on the human phosphorylation S/T predictions reproduces Supplementary Table 16 exactly:
- TT +0.0311 (+0.0273, +0.0351);
- TU −0.0266 (−0.0319, −0.0215);
- UT +0.0236;
- UU +0.0308.

The same run gives the real-data control's matched result (+0.029, row 2), and decision row A for the human task.

## Next

After the END mail for array 4437958: `bash scripts/fetch_mouse.sh` (one Duo push), then
`python -X utf8 -E -P scripts/analyse_mouse.py fetch results`.

## 2026-11-03: results

**Fetch.** `scripts/fetch_mouse.sh`, one connection. Every job finished COMPLETED with exit 0:
- prep 2 min;
- esm 20 min on an A10: 7,069 × 1,280 vectors, mean random-pair cosine 0.891;
- train array 3–16 min per cell, on A100 and A10 GPUs, with no M10.

All 27 json and 54 pred.tsv files are present. The analysis TSV agrees with the per-run JSON AUROCs to within 1e-4.

**Training partitions**, mean over splits:

| | Pure-positive sites | Positive rate |
|---|---|---|
| T | 14.7% | 0.217 |
| Matched | 0.3% | 0.217 |
| U | 0.1% | 0.167 |

**Outcomes** (`scripts/analyse_mouse.py fetch results` → `results/summary.txt`, `results/mouse_results.tsv`). ESM-2 ΔAUROC with 95% joint protein bootstrap intervals:

| Cell | ΔAUROC (95% CI) | Splits > 0 | ΔAUPRC | Permuted ΔAUROC | Real > permuted |
|---|---|---|---|---|---|
| TT | +0.0096 (−0.0172, +0.0332) | 2/3 | −0.0018 | −0.1276 | 3/3 |
| TU | −0.2117 (−0.2245, −0.1982) | 0/3 | −0.4115 | −0.1179 | 0/3 |
| UT | +0.0278 (+0.0218, +0.0335) | 3/3 | | | |
| UU | +0.0206 (+0.0165, +0.0248) | 3/3 | +0.0413 | | |
| MT | +0.0249 (+0.0182, +0.0317) | 3/3 | | | |
| MU | +0.0204 (+0.0163, +0.0254) | 3/3 | | | |

**Pre-specified reading:**
- Decision row **B**: TU is below zero and UU above zero; A does not hold, because the TT interval includes zero.
- Matched control: row **2**. MU is positive in 3/3 splits and close to UU, far from TU.

**Row-B label.** PREREG_mouse.md fixed the label of each row as well as its condition, and it required the label of the observed row to be reported. The label is reported as fixed: "Training-construction effect replicated; the threshold-sampled test does not hide it."
- The first clause describes the result.
- The second clause, read after the results were known, can be taken to say that the threshold-sampled test reveals the harm. It does not: the TT interval includes zero, so TT shows neither a clear benefit nor a harm.
- That clause is recorded as a **post-result interpretive wording correction**. The result is described as a replication of the training-construction effect on an independent source, in which a benefit on threshold-sampled tests was not established. It is not described as a full replication.
- The condition, the observed row and every analysis are unchanged, and the text of PREREG_mouse.md is unchanged.

See LABEL_CORRECTION_rowB_20261103.md and Supplementary Note 9 (Supplementary v12).
