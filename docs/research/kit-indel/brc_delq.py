"""Which base quality does bam-readcount -b apply to a deletion? For reads with a deletion starting
at the KIT site, count how many pass q >= 20 under candidate definitions (anchor base before the
deletion, first base after it, min/mean of both). Compare with bam-readcount -b 20 counts."""
import sys
from collections import Counter

import pysam

S0 = 54723602  # 0-based first deleted base
c = Counter()
for r in pysam.AlignmentFile(sys.argv[1]).fetch("chr4", S0 - 1, S0 + 1):
    if r.mapping_quality < 20 or r.cigartuples is None:
        continue
    rp, qp = r.reference_start, 0
    for op, ln in r.cigartuples:
        if op == 2 and rp == S0:
            q = r.query_qualities
            prev, nxt = q[qp - 1], q[qp]
            k = f"-{ln}"
            c[(k, "n")] += 1
            c[(k, "prev")] += prev >= 20
            c[(k, "next")] += nxt >= 20
            c[(k, "min")] += min(prev, nxt) >= 20
            c[(k, "mean")] += (prev + nxt) / 2 >= 20
            break
        if op in (0, 7, 8):
            rp += ln; qp += ln
        elif op in (2, 3):
            rp += ln
        elif op in (1, 4):
            qp += ln
for k in ("-1", "-2", "-3"):
    print(k, {m: c[(k, m)] for m in ("n", "prev", "next", "min", "mean")})
