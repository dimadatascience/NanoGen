"""Acceptance 3 review: cells with hard flips at >= 2 independent loci between S0 and S6.

Prints each fault cell's per-target (k, m, o) at S0 and S6 and its call at every rung, with the
barcode replaced by F1..Fn. Reads <out>/private/fault_cells.csv (HPC only); print to the job log.
"""
import sys

import pandas as pd

f = pd.read_csv(sys.argv[1] + "/private/fault_cells.csv")
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
f["cell"] = f.cell.map({c: f"F{i + 1}" for i, c in enumerate(f.cell.unique())})
calls = [c for c in f.columns if c.startswith("call_")]
cov = (f[["WT", "MUT", "OTH", "WT_6", "MUT_6", "OTH_6"]].sum(axis=1) > 0)
print(f[cov][["cell", "tid", "locus", "WT", "MUT", "OTH", "WT_6", "MUT_6", "OTH_6", "hard"] + calls]
      .to_string(index=False))
