#!/usr/bin/env python
"""UMI grouping per cell x target bundle: the Levenshtein hybrid merge.

Only primary mapped records (no flag in 0x904) with MAPQ >= min_mapq are used.
A read joins the bundle of every target it reaches: at an SNV, its alignment
covers the site; at an indel, its reference span, deletions included, overlaps
the target +-INDEL_K bp. There is no position term.

Within a bundle, UMI a absorbs b if LD(a, b) <= max_ld and either
  n_a >= ratio * n_b - 1                (directional ratio), or
  n_b <= small_reads and n_a >= n_b     (small clause, kept downhill)
where n is the number of reads per raw UB in the bundle. Clusters grow from the
most abundant UMI through its absorbed UMIs; ties break by count, then by UMI
sequence. The allele is never used. Measured on PALM14E in
docs/research/umi-hybrid-rule.md.

As a script, writes one BAM per cell for bam-readcount -p: a record per read
and bundle, in one read group per cluster, named <target_id>_<cluster ID> and
also set as the MI tag.
"""

import argparse
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import NamedTuple

import pysam

from genotyping import target_id

# Flags of unmapped (0x4), secondary (0x100) and supplementary (0x800) records
NOT_PRIMARY_MAPPED = 0x904
# k: indel bundles reach k bp either side of the target, the same k as the
# tolerant indel ALT match
INDEL_K = 1


def levenshtein(a, b, max_d):
    """Levenshtein distance between a and b, or max_d + 1 if it exceeds max_d.

    Banded DP: only cells within max_d of the diagonal are filled, and the scan
    stops as soon as a whole row exceeds max_d.
    """
    if a == b:
        return 0
    if abs(len(a) - len(b)) > max_d:
        return max_d + 1
    too_far = max_d + 1
    # prev[j] is the distance between a[:i-1] and b[:j]; cells outside the band are too_far
    prev = [j if j <= max_d else too_far for j in range(len(b) + 1)]
    for i in range(1, len(a) + 1):
        lo, hi = max(1, i - max_d), min(len(b), i + max_d)
        cur = [too_far] * (len(b) + 1)
        if i <= max_d:
            cur[0] = i
        for j in range(lo, hi + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1]), too_far)
        if min(cur[max(0, lo - 1) : hi + 1]) > max_d:
            return too_far
        prev = cur
    return prev[-1]


def close_pairs(umis, max_d):
    """Pairs (i, j, d) of distinct UMIs in umis with Levenshtein distance d <= max_d.

    Candidates come from deletion neighbourhoods: two strings within max_d edits
    share a variant obtained by deleting at most max_d bases from each.
    """
    index = defaultdict(set)
    for i, u in enumerate(umis):
        index[u].add(i)
        for n in range(1, max_d + 1):
            for pos in combinations(range(len(u)), n):
                index["".join(c for m, c in enumerate(u) if m not in pos)].add(i)
    candidates = set()
    for members in index.values():
        candidates.update(combinations(sorted(members), 2))
    out = []
    for i, j in sorted(candidates):
        d = levenshtein(umis[i], umis[j], max_d)
        if d <= max_d:
            out.append((i, j, d))
    return out


def absorbs(n_a, n_b, ratio, small_reads):
    """Whether a UMI with n_a reads absorbs a neighbour with n_b reads."""
    return n_a >= ratio * n_b - 1 or (n_b <= small_reads and n_a >= n_b)


def cluster_umis(counts, *, max_ld, ratio, small_reads):
    """Group the raw UBs of one bundle with the hybrid merge.

    counts: raw UB -> reads carrying it (positive). Returns raw UB -> cluster ID;
    IDs are 0, 1, ... in seed order (by count descending, then UB).
    """
    if any(n <= 0 for n in counts.values()):
        raise ValueError("Reads per raw UB must be positive")
    umis = sorted(counts, key=lambda u: (-counts[u], u))
    neighbours = defaultdict(list)
    for i, j, _ in close_pairs(umis, max_ld):
        neighbours[i].append(j)
        neighbours[j].append(i)

    cluster = {}
    n_clusters = 0
    for seed in range(len(umis)):
        if seed in cluster:
            continue
        cluster[seed] = n_clusters
        stack = [seed]
        while stack:
            a = stack.pop()
            for b in sorted(neighbours[a]):
                if b not in cluster and absorbs(counts[umis[a]], counts[umis[b]], ratio, small_reads):
                    cluster[b] = n_clusters
                    stack.append(b)
        n_clusters += 1
    return {umis[i]: c for i, c in cluster.items()}


##


# Bundles


class Target(NamedTuple):
    """A target and its bundle: the 0-based, half-open reference interval
    [lo, hi) that a read's reference span must overlap."""

    target_id: str
    chrom: str
    lo: int
    hi: int

    @classmethod
    def from_bed(cls, chrom, start, ref, alt, k=INDEL_K):
        """start is the BED start, the 1-based bam-readcount position: the
        SNV site, the first deleted base or the insertion anchor."""
        pos = int(start) - 1
        if alt.startswith("-"):
            lo, hi = pos - k, pos + len(alt) - 1 + k
        elif alt.startswith("+"):
            lo, hi = pos - k, pos + 1 + k
        else:
            lo, hi = pos, pos + 1
        return cls(target_id(chrom, start, ref, alt), chrom, lo, hi)


def read_targets(bed):
    """Targets of a BED (chr, start, end, ref, alt, gene), one per target ID, in file order."""
    targets = {}
    with open(bed) as f:
        for line in f:
            if line.strip():
                chrom, start, _, ref, alt = line.rstrip("\n").split("\t")[:5]
                t = Target.from_bed(chrom, start, ref, alt)
                targets.setdefault(t.target_id, t)
    return list(targets.values())


def is_primary_mapped(read, min_mapq):
    """Whether a record enters grouping: primary, mapped and MAPQ >= min_mapq."""
    return not read.flag & NOT_PRIMARY_MAPPED and read.mapping_quality >= min_mapq


def in_bundle(read, target):
    """Whether a read's reference span, deletions included, overlaps the target's bundle."""
    return (read.reference_name == target.chrom
            and read.reference_start < target.hi
            and read.reference_end > target.lo)


def group_reads(reads, targets, *, min_mapq, max_ld, ratio, small_reads):
    """
    Group one cell's reads per target with the hybrid merge.

    Records that are not primary mapped, below min_mapq or without a UB tag are
    dropped. A read in several bundles is grouped in each, on that bundle's
    reads per UB.

    Yields:
        Tuple[str, int, pysam.AlignedSegment]: (target ID, cluster ID, read) per
        read and bundle, target by target in the order given; cluster IDs
        restart at 0 per target (see cluster_umis).
    """
    bundles = {t.target_id: [] for t in targets}
    for read in reads:
        if not is_primary_mapped(read, min_mapq) or not read.has_tag("UB"):
            continue
        for t in targets:
            if in_bundle(read, t):
                bundles[t.target_id].append(read)
    for tid, members in bundles.items():
        counts = Counter(r.get_tag("UB") for r in members)
        clusters = cluster_umis(counts, max_ld=max_ld, ratio=ratio, small_reads=small_reads)
        for r in members:
            yield tid, clusters[r.get_tag("UB")], r


def write_grouped_bam(input, bed, output, *, min_mapq, max_ld, ratio, small_reads):
    """
    Write a cell's grouped reads for bam-readcount -p: one record per read and
    bundle, in a read group per cluster named <target_id>_<cluster ID> (ID and
    LB), also set as the MI tag. The output is coordinate-sorted and indexed.
    """
    with pysam.AlignmentFile(str(input), "rb", check_sq=False) as f:
        header = f.header.to_dict()
        grouped = list(group_reads(f, read_targets(bed), min_mapq=min_mapq, max_ld=max_ld,
                                   ratio=ratio, small_reads=small_reads))

    # Read groups by target, then cluster ID
    order = {tid: i for i, tid in enumerate(dict.fromkeys(tid for tid, _, _ in grouped))}
    clusters = sorted({(order[tid], c, tid) for tid, c, _ in grouped})
    cell = Path(input).name.split(".")[0]
    header["HD"] = {**header.get("HD", {"VN": "1.6"}), "SO": "coordinate"}
    header["RG"] = [{"ID": f"{tid}_{c}", "LB": f"{tid}_{c}", "SM": cell} for _, c, tid in clusters]
    out_header = pysam.AlignmentHeader.from_dict(header)

    records = []
    for tid, c, read in grouped:
        record = pysam.AlignedSegment.from_dict(read.to_dict(), out_header)
        record.set_tag("RG", f"{tid}_{c}", "Z")
        record.set_tag("MI", f"{tid}_{c}", "Z")
        records.append(record)
    records.sort(key=lambda r: (r.reference_id, r.reference_start))

    with pysam.AlignmentFile(str(output), "wb", header=out_header) as out:
        for record in records:
            out.write(record)
    pysam.index(str(output))


##


def _parse_args():
    parser = argparse.ArgumentParser(description="Group one cell's reads per target by UMI.")
    parser.add_argument("-i", "--input", required=True, help="Per-cell BAM with UB tags.")
    parser.add_argument("-b", "--bed", required=True,
                        help="Targets: chr, start, end, ref, alt, gene (tab-separated).")
    parser.add_argument("-o", "--output", required=True, help="Grouped, indexed BAM to write.")
    parser.add_argument("--min_mapq", type=int, default=20, help="Minimum MAPQ.")
    parser.add_argument("--max_ld", type=int, default=2, help="Maximum UMI Levenshtein distance.")
    parser.add_argument("--ratio", type=float, default=2,
                        help="Directional ratio: a absorbs b if n_a >= ratio * n_b - 1.")
    parser.add_argument("--small_reads", type=int, default=3,
                        help="Small clause: a absorbs b if n_b <= small_reads and n_a >= n_b.")
    return parser.parse_args()


def main():
    args = _parse_args()
    write_grouped_bam(args.input, args.bed, args.output, min_mapq=args.min_mapq,
                      max_ld=args.max_ld, ratio=args.ratio, small_reads=args.small_reads)


if __name__ == "__main__":
    main()
