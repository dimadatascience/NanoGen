"""What are the giant UMI groups at the indel target site? Size distribution of bam-readcount libraries at the
indel target, and for groups with >= 50 reads: UMI composition (missing / low complexity / dominant
base), cells involved, and new-rule labels. Aggregates only."""
import glob
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from kit_indel import NMIN_MAIN, main_label, new_counts, parse_line, uniform_label  # noqa: E402

cells_dir = sys.argv[1]
chrom, pos, ref, alt = next((f[0], int(f[1]), f[3], f[4]) for f in (l.rstrip("\n").split("\t") for l in open(sys.argv[2]))
                            if f[4][0] in "+-")
rows = []
for cd in sorted(glob.glob(os.path.join(cells_dir, "*"))):
    cp = os.path.join(cd, "allcounts.count")
    if not os.path.exists(cp):
        continue
    with open(cp) as fh:
        for line in fh:
            f = line.split("\t", 3)
            if f[0] == chrom and int(f[1]) == pos:
                for lb, e in parse_line(line).items():
                    umi = lb.rsplit("_", 2)[0]
                    c = new_counts(e, ref, alt, 1)
                    rows.append({"cell": os.path.basename(cd), "umi": umi, "n": sum(c.values()),
                                 "alt_frac": c["ALT"] / max(sum(c.values()), 1),
                                 "main": main_label(e, ref, alt),
                                 "new": uniform_label(c, NMIN_MAIN, NMIN_MAIN, NMIN_MAIN)})
D = pd.DataFrame(rows)
print("libraries", len(D), "reads", D.n.sum())
bins = [0, 2, 5, 10, 50, 100, 1000, 10**7]
D["size"] = pd.cut(D.n, bins, right=False)
print(D.groupby("size", observed=True).agg(libs=("n", "size"), reads=("n", "sum"),
                                           alt_frac_median=("alt_frac", "median")).to_string())


def kind(u):
    if u in ("None", "", "nan"):
        return "missing"
    top = Counter(u).most_common(1)[0][1] / len(u)
    return "low_complexity(top base >= 75%)" if top >= .75 else ("dominant base 50-75%" if top >= .5 else "normal")


B = D[D.n >= 50].copy()
B["umi_kind"] = B.umi.map(kind)
print("\ngiant (>=50 reads):", len(B), "libraries,", B.cell.nunique(), "cells,",
      f"{B.n.sum() / D.n.sum():.1%} of site reads; umi length", Counter(B.umi.str.len()).most_common(3))
print(B.groupby("umi_kind").agg(libs=("n", "size"), reads=("n", "sum")).to_string())
print("top base of low-complexity UMIs:", Counter(Counter(u).most_common(1)[0][0] for u in B.umi if kind(u).startswith("low")))
print("labels main -> new:", B.groupby(["main", "new"]).size().to_dict())
print("giant libraries per cell (cells with any):", B.groupby("cell").size().describe().round(1).to_dict())
print("normal-size (<50) UMIs low complexity frac:", np.mean([kind(u).startswith("low") for u in D[D.n < 50].umi]))
