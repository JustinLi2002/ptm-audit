# Analysed data

Two releases from http://deepmvp.ptmax.org were analysed.

## Individually posted .tsv files

Downloaded [DATE]. File timestamps 2025-01-11. 31-residue windows.
Confirmed by the authors (personal communication, July 2026) to predate the
final pipeline: methylation K and R are merged in these files, and the test
partitions are class-balanced rather than split as described in the paper.

    md5sum retrain/*.tsv > checksums_tsv.txt

## all_data.tar.gz

Downloaded [DATE]. 61-residue windows, protein-level split, natural class
ratios, methylation K and R separate. Test partitions filtered at 70/80/90%
peptide identity to the training set.

    md5sum allbig/all_data/*.tsv > checksums_alldata.txt
