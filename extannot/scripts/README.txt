External-annotation check (PREREG.md). Run from any directory; the data folders are passed as arguments.

1. python -X utf8 -E -P scripts/build_sites.py data results_DIR
       parses UniProt (data/uniprot), CPLM (data/cplm) and, when present, PhosphoSitePlus (data/psp)
2. python -X utf8 -E -P scripts/extannot_check.py data/zenodo_x results_DIR/sites.tsv results_DIR
       Step 1 counts, Step 2 rescoring, enrichment, decision rule
3. python -X utf8 -E -P scripts/make_results_md.py results_DIR > tables.md

data/zenodo_x holds the four archives of Zenodo record 21670043, extracted; checksums in zenodo_md5.txt.
fetch_uniprot.py fetched data/uniprot/swissprot_human_ptm_2026_03.tsv.
