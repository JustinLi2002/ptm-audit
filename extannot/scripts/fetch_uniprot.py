"""Page through UniProt search (the stream endpoint returned HTTP 500 on 2026-11-02) and write one TSV:
human reference proteome, reviewed entries, with modified-residue, glycosylation and cross-link features."""
import re
import sys
import urllib.request

OUT = sys.argv[1]
url = ("https://rest.uniprot.org/uniprotkb/search?format=tsv&size=500"
       "&fields=accession,ft_mod_res,ft_carbohyd,ft_crosslnk"
       "&query=%28proteome%3AUP000005640%29%20AND%20%28reviewed%3Atrue%29")
rows, header, release, total = 0, None, None, None
with open(OUT, "w", encoding="utf-8", newline="\n") as out:
    while url:
        with urllib.request.urlopen(url, timeout=120) as r:
            release = release or r.headers.get("X-UniProt-Release")
            total = total or int(r.headers.get("X-Total-Results"))
            lines = r.read().decode("utf-8").splitlines()
            m = re.search(r'<([^>]+)>;\s*rel="next"', r.headers.get("Link") or "")
            url = m.group(1) if m else None
        if header is None:
            header = lines[0]
            out.write(header + "\n")
        else:
            assert lines[0] == header
        for l in lines[1:]:
            out.write(l + "\n")
            rows += 1
assert rows == total, (rows, total)
print(f"release {release}: {rows} entries written to {OUT}")
