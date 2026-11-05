"""Figure 1 (v7). As fig1_construction.py, with the design note naming the three protein-disjoint
train-validation-test splits (v7 terminology: a split is one of the three replicates; a partition is its
training, validation or test part). Runs the original script with two strings replaced; that file is not changed."""
from pathlib import Path

src = (Path(__file__).resolve().parent / "fig1_construction.py").read_text(encoding="utf-8")
EDITS = [('"Same three protein-disjoint partitions in every cell;\\n"\n'
          '                               "each trained ensemble is scored on both test constructions."',
          '"The same three protein-disjoint train–validation–test\\n"\n'
          '                               "splits in every cell; each trained ensemble\\n"\n'
          '                               "is scored on both test constructions."'),
         ('save(fig, "Figure1_construction_coverage")', 'save(fig, "Figure1_construction_coverage_v7")')]
for old, new in EDITS:
    assert src.count(old) == 1, old
    src = src.replace(old, new)
exec(compile(src, "fig1_construction_v7", "exec"), {"__name__": "__main__", "__file__": __file__})
