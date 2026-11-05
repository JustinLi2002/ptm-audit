External-annotation check (Supplementary Note 6, Supplementary Table 15, Methods).

Asks whether test negatives recorded as sites in external databases explain the harm on the unrestricted test. Uses the
deposited predictions only; nothing is retrained. Design and decision rule: PREREG.md. Step-by-step commands:
scripts/README.txt.

results/               the final run: UniProt, CPLM and PhosphoSitePlus phosphorylation
results_uniprot_cplm/  the earlier run without PhosphoSitePlus

Only aggregate tables are here. The two site-level lists the scripts write, sites.tsv and flagged_test_negatives.tsv,
derive from PhosphoSitePlus and CPLM and are not redistributed; neither are the database files themselves.
scripts/build_sites.py and scripts/list_flagged.py regenerate both lists from the files below.

Obtaining the databases

UniProtKB/Swiss-Prot, release 2026_03
    scripts/fetch_uniprot.py repeats the query. It returns the current release, so the result may differ from the
    release used here.

CPLM 4.0 (Zhang et al., Nucleic Acids Res 2022;50:D451-9)
    The Homo sapiens archive from the CPLM download page, site last updated 10 Feb 2025.

PhosphoSitePlus (Hornbeck et al.)
    Phosphorylation_site_dataset, a copy downloaded on 2025-06-06 (MD5 585d8cecaa587b80f545d68d136db617). Current
    PhosphoSitePlus terms require registration and may not allow downloading without an agreement with Cell Signaling
    Technology, so readers need to obtain access themselves. The check can be run without it, as in
    results_uniprot_cplm/.

data/zenodo_x in the scripts holds the extracted archives of the data record (first version doi:10.5281/zenodo.21670043);
checksums are in scripts/zenodo_md5.txt.
