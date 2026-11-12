#!/usr/bin/env python3
"""Exit 0 iff two summarise tables match: same header/rows, non-numeric fields equal, numbers within 1e-9
(the last float bit varies with summation order)."""
import sys, math, csv
a = list(csv.reader(open(sys.argv[1]), delimiter="\t")); b = list(csv.reader(open(sys.argv[2]), delimiter="\t"))
ok = len(a) == len(b) and a[0] == b[0]
worst = 0.0
for ra, rb in zip(a[1:], b[1:]):
    for x, y in zip(ra, rb):
        if x == y: continue
        try:
            fx, fy = float(x), float(y)
        except ValueError:
            ok = False; continue
        if math.isnan(fx) and math.isnan(fy): continue
        d = abs(fx - fy); worst = max(worst, d)
        if not d <= 1e-9: ok = False
print(f"match={ok} max_abs_diff={worst:.2e}")
sys.exit(0 if ok else 1)
