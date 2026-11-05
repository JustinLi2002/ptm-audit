"""Externally annotated human modification sites, one row per (source, tier, accession, position, residue, kind),
following PREREG.md (ptm-audit/extannot_20261102). Reads UniProt (TSV from fetch_uniprot.py), CPLM (Homo sapiens
archive) and, when present, the PhosphoSitePlus site datasets; writes sites.tsv and a parse report. CPLM rows carry
the 15-residue window from CPLM's own sequence, to check sequence drift against our reference.

kind is the modification class used by the task mapping: phospho, acetyl, methyl_K, methyl_R, ubiquitin, sumo, nglyc.

Usage: python -I build_sites.py DATA_DIR OUT_DIR"""
import csv
import gzip
import io
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
rows, report = [], []

# ---------------------------------------------------------------- UniProt
EXPERIMENTAL = ("ECO:0000269", "ECO:0007744")
SIMILARITY = "ECO:0000250"
FEAT = re.compile(r'(MOD_RES|CARBOHYD|CROSSLNK) (\d+)(?:\.\.(\d+))?;\s*/note="([^"]*)"(?:;\s*/evidence="([^"]*)")?')


def uniprot_kind(ftype, note):
    n = note.split(";")[0].strip()
    if ftype == "MOD_RES":
        if n == "Phosphoserine":
            return "phospho", "S"
        if n == "Phosphothreonine":
            return "phospho", "T"
        if n == "Phosphotyrosine":
            return "phospho", "Y"
        if n == "N6-acetyllysine":
            return "acetyl", "K"
        if re.fullmatch(r"N6(,N6){0,2}-(tri|di)?methyllysine", n):
            return "methyl_K", "K"
        if re.fullmatch(r"(Omega-N-methylarginine|Asymmetric dimethylarginine|Symmetric dimethylarginine|"
                        r"Dimethylated arginine|N5-methylarginine)", n):
            return "methyl_R", "R"
        return None
    if ftype == "CARBOHYD":
        if n.startswith("N-linked") and n.endswith("asparagine"):
            return "nglyc", "N"
        return None
    if ftype == "CROSSLNK" and n.startswith("Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in"):
        partner = n.rsplit(" in ", 1)[1].rstrip(")")
        if partner == "ubiquitin":
            return "ubiquitin", "K"
        if re.fullmatch(r"SUMO[123]?", partner):
            return "sumo", "K"
    return None


up = DATA / "uniprot" / "swissprot_human_ptm_2026_03.tsv"
seen_notes = Counter()
with open(up, encoding="utf-8") as fh:
    rd = csv.reader(fh, delimiter="\t")
    header = next(rd)
    assert header == ["Entry", "Modified residue", "Glycosylation", "Cross-link"], header
    for acc, *cols in rd:
        for col in cols:
            for ftype, a, b, note, ev in FEAT.findall(col):
                if b:  # ranges are intrachain links, never a single modified residue
                    continue
                k = uniprot_kind(ftype, note)
                if k is None:
                    seen_notes[(ftype, note.split(";")[0])] += 1
                    continue
                if any(c in ev for c in EXPERIMENTAL):
                    tier = "primary"
                elif SIMILARITY in ev:
                    tier = "similarity"
                else:
                    continue
                rows.append(("UniProt", tier, acc, int(a), k[1], k[0], ""))
report.append(f"UniProt {up.name}: {sum(r[0] == 'UniProt' for r in rows)} sites kept")
report += [f"  unmapped note (top): {f} {n!r} x{c}" for (f, n), c in seen_notes.most_common(25)]

# ---------------------------------------------------------------- CPLM (lysine only)
CPLM_KIND = {"Acetylation": "acetyl", "Methylation": "methyl_K", "Ubiquitination": "ubiquitin",
             "Sumoylation": "sumo"}
bad = Counter()
with zipfile.ZipFile(DATA / "cplm" / "Homo_sapiens.zip") as z:
    (name,) = z.namelist()
    with z.open(name) as fh:
        for line in io.TextIOWrapper(fh, encoding="utf-8"):
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or f[3] not in CPLM_KIND:
                continue
            assert f[5] == "Homo sapiens", f[5]
            pos, seq = int(f[2]), f[6]
            if not (0 < pos <= len(seq)) or seq[pos - 1] != "K":
                bad[f[3]] += 1
                continue
            win = seq[pos - 8:pos + 7] if pos >= 8 and pos + 7 <= len(seq) else ""
            rows.append(("CPLM", "primary", f[1], pos, "K", CPLM_KIND[f[3]], win))
report.append(f"CPLM {name}: {sum(r[0] == 'CPLM' for r in rows)} sites kept; not K in its own sequence: {dict(bad)}")

# ---------------------------------------------------------------- PhosphoSitePlus (optional until downloaded)
PSP_FILES = {"Phosphorylation_site_dataset": "phospho", "Acetylation_site_dataset": "acetyl",
             "Methylation_site_dataset": "methyl", "Ubiquitination_site_dataset": "ubiquitin",
             "Sumoylation_site_dataset": "sumo"}
psp_dir = DATA / "psp"
for stem, kind in PSP_FILES.items():
    path = next((p for p in sorted(psp_dir.glob(stem + "*")) if p.is_file()), None) if psp_dir.exists() else None
    if path is None:
        report.append(f"PhosphoSitePlus {stem}: not present")
        continue
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    # the header follows a licence preamble in some releases, and its column order varies
    h = next(i for i, l in enumerate(lines) if "\tACC_ID\t" in l and "\tMOD_RSD\t" in l)
    rd = csv.DictReader(lines[h:], delimiter="\t")
    n0 = sum(r[0] == "PhosphoSitePlus" for r in rows)
    skipped = Counter()
    for r in rd:
        if r["ORGANISM"] != "human":
            continue
        acc = r["ACC_ID"]
        if "-" in acc:  # isoform-specific accession
            skipped["isoform"] += 1
            continue
        m = re.fullmatch(r"([A-Z])(\d+)-\w+", r["MOD_RSD"].strip())
        if not m:
            skipped["unparsed MOD_RSD"] += 1
            continue
        res, pos = m.group(1), int(m.group(2))
        k = kind
        if kind == "methyl":
            k = {"K": "methyl_K", "R": "methyl_R"}.get(res)
            if k is None:
                skipped["methyl residue not K/R"] += 1
                continue
        lt = int(r["LT_LIT"] or 0) if (r.get("LT_LIT") or "").strip().isdigit() else 0
        rows.append(("PhosphoSitePlus", "primary", acc, pos, res, k, ""))
        if lt >= 1:
            rows.append(("PhosphoSitePlus", "low_throughput", acc, pos, res, k, ""))
    report.append(f"PhosphoSitePlus {path.name}: {sum(r[0] == 'PhosphoSitePlus' for r in rows) - n0} rows kept "
                  f"(primary + low-throughput); skipped {dict(skipped)}")

rows = sorted(set(rows))
with open(OUT / "sites.tsv", "w", encoding="utf-8", newline="\n") as fh:
    fh.write("source\ttier\taccession\tposition\tresidue\tkind\twindow15\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
summary = Counter((r[0], r[1], r[5]) for r in rows)
report.append("unique sites by source, tier, kind:")
report += [f"  {s}\t{t}\t{k}\t{n}" for (s, t, k), n in sorted(summary.items())]
(OUT / "sites_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
print("\n".join(report))
