"""Merge gate (#41): diagnose clusters whose usable-read counts differ between v2 and the emulation.

For each cluster in <gate>/private/cluster_disagreements.csv.gz, re-groups the cell with v2's
umi_grouping, writes each read of the cluster to its own BAM, runs bam-readcount (-b 0 at indels,
-b 20 at SNVs) on it and compares with the emulation's per-read allele. Prints the reads of those cells' bundles that differ
as CIGAR operations within 5 bp of the site (no names, sequences or barcodes).
Usage: python diagnose_counts.py <gate dir> <cell BAM dir> <v2 bin> <bed> <ref> [--with-others]
"""
import os
import subprocess
import sys
import tempfile

import pandas as pd
import pysam

gate, cells, v2_bin, bed, ref = sys.argv[1:6]
# --with-others: keep the other targets' records (same reads, other bundles) in the BAM, as in v2
WITH_OTHERS = "--with-others" in sys.argv
sys.path.insert(0, v2_bin)
import count_consensus as cc  # noqa: E402
import umi_grouping as ug  # noqa: E402
from genotyping import target_id  # noqa: E402

from kit_indel import read_info  # noqa: E402
from measure_umi import load_targets  # noqa: E402
from merge_gate import BQ, indel_allele  # noqa: E402

OPS = "MIDNSHP=X"
T = load_targets(bed)
D = pd.read_csv(os.path.join(gate, "private", "cluster_disagreements.csv.gz"))
D = D[D.why == "counts"].drop_duplicates(["cell", "tid"])
fasta = pysam.FastaFile(ref)
for n, row in enumerate(D.itertuples(index=False)):
    t = T[T.tid == row.tid].iloc[0]
    tid = target_id(t.chr, t.start, t.ref, t.alt)
    tmp = tempfile.mkdtemp()
    grouped = os.path.join(tmp, "g.bam")
    ug.write_grouped_bam(os.path.join(cells, f"{row.cell}.bam"), bed, grouped, min_mapq=20, max_ld=2,
                         ratio=2, small_reads=3)
    w0 = t.start - 1 - 50
    seq = fasta.fetch(t.chr, w0, t.start - 1 + 50).upper()
    with pysam.AlignmentFile(grouped) as f:
        # every read of the cell's bundle at this target (a cell may list several clusters)
        reads = [r for r in f.fetch(t.chr, t.start - 2, t.start + len(t.alt))
                 if ug.read_group_target(r.get_tag("MI")) == tid]
        others = [r for r in f.fetch() if ug.read_group_target(r.get_tag("MI")) != tid] if WITH_OTHERS else []
        header = f.header.to_dict()
        header["RG"] = header["RG"] + [{"ID": f"r{i}", "LB": f"r{i}", "SM": "x"} for i in range(len(reads))]
        one = os.path.join(tmp, "one.bam")
        recs = []
        out_h = pysam.AlignmentHeader.from_dict(header)
        for i, r in enumerate(reads):
            rec = pysam.AlignedSegment.from_dict(r.to_dict(), out_h)
            rec.set_tag("RG", f"r{i}", "Z")
            recs.append(rec)
        recs += [pysam.AlignedSegment.from_dict(r.to_dict(), out_h) for r in others]
        recs.sort(key=lambda r: (r.reference_id, r.reference_start))
        with pysam.AlignmentFile(one, "wb", header=out_h) as o:
            for rec in recs:
                o.write(rec)
        pysam.index(one)
        site = os.path.join(tmp, "site.bed")
        with open(site, "w") as s:
            s.write(f"{t.chr}\t{t.start}\t{t.end}\t{t.ref}\t{t.alt}\t{t.gene}\n")
        line = subprocess.run(["bam-readcount", "-f", ref, "--site-list", site, "-b",
                               "0" if t.kind == "indel" else str(BQ), one, "-p"],
                              capture_output=True, text=True, check=True).stdout.splitlines()
        per_read = {}
        if line:
            for lib, entries in cc.parse_readcount_line(line[0])[3]:
                if not lib.startswith("r"):
                    continue
                per_read[int(lib[1:])] = cc.allele_counts(entries, t.ref, t.alt)
        diffs = 0
        for i, r in enumerate(reads):
            v2 = per_read.get(i, (0, 0, 0))
            raw = read_info(r, t, seq, w0)[0]
            emu = indel_allele(raw, t)
            emu = tuple(int(emu == x) for x in ("REF", "ALT", "OTH"))
            if v2 != emu:
                diffs += 1
                # CIGAR ops overlapping site +-5 bp
                ops, pos = [], r.reference_start
                for op, ln in r.cigartuples:
                    end = pos + (ln if op in (0, 2, 3, 7, 8) else 0)
                    if end >= t.start - 6 and pos <= t.start + 5:
                        ops.append(f"{ln}{OPS[op]}@{pos - (t.start - 1):+d}")
                    pos = end
                print(f"cluster {n + 1} ({row.tid}, {row.reads} reads): read v2 r/a/o {v2} emu {emu} "
                      f"site entries {raw} cigar near site {' '.join(ops)}")
        print(f"cluster {n + 1}: {diffs} reads differ")
