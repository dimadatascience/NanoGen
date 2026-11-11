#!/usr/bin/env python3
"""Synthetic smoke test of power_check.py: a tiny BAM with known UMI clusters.

usage: smoke.py OUTDIR   (GENOTYPING_DIR and PYTHONPATH set as for power_check.py)
"""
import os
import subprocess
import sys

import pysam

d = sys.argv[1]
os.makedirs(d, exist_ok=True)
here = os.path.dirname(os.path.abspath(__file__))
C1, C2, C3 = "A" * 16, "C" * 16, "G" * 16
# SNV at 1-based 101 (0-based 100); deletion of 2 bases at 1-based 201 (0-based 200..201)
open(f"{d}/t.bed", "w").write("chr1\t101\t101\tC\tT\tG1\nchr1\t201\t201\tAG\t-AG\tG2\n")
hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
reads = []


def add(cb, ub, start, cigar, q=30, flag=0, mapq=60):
    reads.append((start, cb, ub, cigar, q, flag, mapq))


# C1 / SNV: UMI U1 x3 good (cluster 1 -> 1); U2 x3 with one low-quality read (2 usable -> 0);
#           U1-like U1x (LD 1, one read) merges into U1 (small clause) -> still 1 cluster counted
for _ in range(3):
    add(C1, "AAAAAAAAAAAA", 50, [(0, 100)])
add(C1, "AAAAAAAAAAAC", 50, [(0, 100)])
add(C1, "TTTTTTTTTTTT", 50, [(0, 100)])
add(C1, "TTTTTTTTTTTT", 50, [(0, 100)])
add(C1, "TTTTTTTTTTTT", 50, [(0, 100)], q=5)
# C2 / SNV: one UMI x3 reads, one secondary (ignored) -> 1 ; plus a read spliced over the site (N) -> not in bundle
for _ in range(3):
    add(C2, "GGGGGGGGGGGG", 60, [(0, 80)])
add(C2, "GGGGGGGGGGGG", 60, [(0, 80)], flag=256)
add(C2, "CCCCCCCCCCCC", 60, [(0, 20), (3, 50), (0, 20)])
# C2 / deletion: UMI with 3 reads carrying a deletion -> 1 ; MAPQ 5 read ignored
for _ in range(3):
    add(C2, "ACACACACACAC", 180, [(0, 20), (2, 2), (0, 30)])
add(C2, "ACACACACACAC", 180, [(0, 20), (2, 2), (0, 30)], mapq=5)
# C3 (not in QC/negatives) -> ignored
for _ in range(3):
    add(C3, "GGGGGGGGGGGG", 60, [(0, 80)])

with pysam.AlignmentFile(f"{d}/t.bam", "wb", header=hdr) as out:
    for i, (start, cb, ub, cig, q, flag, mapq) in enumerate(sorted(reads, key=lambda r: r[0])):
        a = pysam.AlignedSegment()
        a.query_name = f"r{i}"
        a.reference_id, a.reference_start, a.cigartuples = 0, start, cig
        ln = sum(l for op, l in cig if op in (0, 1, 4))
        a.query_sequence = "C" * ln
        a.query_qualities = pysam.qualitystring_to_array(chr(q + 33) * ln)
        a.flag, a.mapping_quality = flag, mapq
        a.set_tag("CB", cb)
        a.set_tag("UB", ub)
        out.write(a)
pysam.index(f"{d}/t.bam")
open(f"{d}/qc.txt", "w").write(f"{C1}\n{C2}\n")
open(f"{d}/neg.csv", "w").write(f"sample,cell\nX,{C1}\n")
open(f"{d}/lab.csv", "w").write(f'"cell16","ct"\n"{C2}","Leukemic"\n')
py = [sys.executable, f"{here}/power_check.py"]
subprocess.run(py + ["coverage", "--bam", f"{d}/t.bam", "--bed", f"{d}/t.bed", "--qc", f"{d}/qc.txt",
                     "--negatives", f"{d}/neg.csv", "--out", f"{d}/cells.tsv"], check=True)
got = {(r.split("\t")[0], r.split("\t")[2]): int(r.split("\t")[3]) for r in open(f"{d}/cells.tsv").read().split("\n")[1:] if r}
want = {("T1", C1): 1, ("T1", C2): 1, ("T2", C2): 1}
assert got == want, (got, want)
subprocess.run(py + ["summarise", "--sample", "X", "--cells", f"{d}/cells.tsv", "--qc", f"{d}/qc.txt",
                     "--negatives", f"{d}/neg.csv", "--labels", f"{d}/lab.csv", "--n-targets", "2",
                     "--out", f"{d}/sum.tsv"], check=True)
sys.path.insert(0, here)
import power_check as P  # noqa: E402
pn = P.power_numbers(30.0)
assert 1.3 < pn["r80"] < 1.6 and pn["fail_r15"] > pn["fail_r1"], pn
assert P.power_numbers(0.0)["r80"] == float("inf")
print("smoke OK", got, pn, f"E*={P.e_star():.1f}")
