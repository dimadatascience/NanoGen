"""Compare the subset re-run against the main@c442a32 baseline (aggregates only)."""
import sys
import pandas as pd

base = pd.read_csv(sys.argv[1], index_col=0)
new = pd.read_csv(sys.argv[2], index_col=0)
key = ["cell", "chr", "start", "alt"]
m = base.merge(new, on=key, how="left", suffixes=("_b", "_n"))
for c in ["WT", "MUT", "MIS"]:
    m[f"{c}_n"] = m[f"{c}_n"].fillna(0)
m["genotype_b"] = m["genotype_b"].fillna("NA")
m["genotype_n"] = m["genotype_n"].fillna("NA")
same = ((m.WT_b == m.WT_n) & (m.MUT_b == m.MUT_n) & (m.MIS_b == m.MIS_n))
print("rows", len(m), "count-identical", int(same.sum()), f"({same.mean():.4%})")
print("UMI totals base:", int(m.WT_b.sum()), int(m.MUT_b.sum()), int(m.MIS_b.sum()),
      " rerun:", int(m.WT_n.sum()), int(m.MUT_n.sum()), int(m.MIS_n.sum()))
print(pd.crosstab(m.genotype_b, m.genotype_n, margins=True))
