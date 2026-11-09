"""Synthetic `main` input for the #67 smoke test: main-format counts (with MIS) from the synthetic v3 table.
The output is then called by main's own bin/long_reads/genotyping.py (see run.sbatch)."""
import sys

import numpy as np
import pandas as pd

rng = np.random.default_rng(2)
v3 = pd.read_csv(sys.argv[1])
n = (v3["WT"] + v3["MUT"]).values
pd.DataFrame(dict(chr="chr1", start=0, end=0, ref="A", alt="T", gene=v3["gene"], WT=v3["WT"], MUT=v3["MUT"],
                  MIS=rng.poisson(0.05 * n + 0.02), barcode=v3["cell"].str.slice(2, 18))).to_csv(sys.argv[2], index=False)
