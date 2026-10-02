"""Tiny synthetic fixture to smoke-test measure_umi.py (no real data)."""
import os, random, pysam
random.seed(0)
os.makedirs("synth/cells", exist_ok=True)
with open("synth/t.bed", "w") as f:
    f.write("chrS\t50001\t50001\tA\tC\tG1\n")
    f.write("chrS\t50005\t50005\tG\t-TTA\tG2\n")
hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrS", "LN": 100000}]}
def mk(f, name, start, ub, bx, ug, base, flag=0, dele=False):
    a = pysam.AlignedSegment(f.header)
    a.query_name = name; a.reference_id = 0; a.reference_start = start; a.mapping_quality = 60; a.flag = flag
    n1 = 50004 - start + 1  # through anchor (0-based 50004)
    if dele:
        a.cigartuples = [(0, n1), (2, 3), (0, 100)]; L = n1 + 100
    else:
        a.cigartuples = [(0, 300)]; L = 300
    s = [random.choice("ACGT") for _ in range(L)]
    s[50000 - start] = base; s[50004 - start] = "G"
    a.query_sequence = "".join(s); a.query_qualities = pysam.qualitystring_to_array("I" * L)
    a.set_tag("UB", ub); a.set_tag("BX", bx); a.set_tag("UG", ug)
    return a
for cell, recs in {
    "CELLAAAA": [("r1", 49950, "AAAACCCCGGGG", "C", 0, False), ("r2", 49960, "AAAACCCCGGGG", "C", 0, False),
                 ("r3", 50000, "AAAACCCCGGGG", "C", 0, False), ("r4", 49955, "AAAACCCCGGGG", "C", 256, False),
                 ("r5", 49970, "AAAACCCCGGGT", "C", 0, False), ("r6", 49950, "TTTTCCCCAAAA", "A", 0, True),
                 ("r7", 49951, "TTTTCCCCAAAA", "A", 0, True), ("r8", 49952, "TTTTCCCCAAA", "A", 0, False)],
    "CELLCCCC": [("s1", 49950, "GGGGTTTTAAAA", "A", 0, False), ("s2", 49950, "GGGGTTTTAAAA", "A", 0, False),
                 ("s3", 49950, "GGGGTTTAAAAC", "A", 0, False)]}.items():
    d = f"synth/cells/{cell}"; os.makedirs(d, exist_ok=True)
    with pysam.AlignmentFile(f"{d}/filtered_input.bam", "wb", header=hdr) as f:
        for i, (n, st, ub, b, fl, de) in enumerate(sorted(recs, key=lambda x: x[1])):
            f.write(mk(f, n, st, ub, ub, i, b, fl, de))
    pysam.index(f"{d}/filtered_input.bam")
print("ok")
