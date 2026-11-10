#!/usr/bin/env python
"""UMI grouping per cell x target bundle: the Levenshtein hybrid merge.

Only primary mapped records (no flag in 0x904) with MAPQ >= min_mapq are used.
A read joins the bundle of every target it reaches: at an SNV, its alignment
covers the site; at an indel, its reference span, deletions included, overlaps
the target +-INDEL_K bp. A skipped intron (N) is not part of either. There is
no position term.

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
# k: indel bundles reach k bp either side of the target, and an indel entry is
# ALT-like if its length is within k bp of the target's (count_consensus.py)
INDEL_K = 1
# PROTOTYPE (position-tolerant indel matching, never merged): W, the bp either
# side of the anchor where a tolerant indel's left end may sit. -1 switches the
# rule off and gives the genotype-v3 counts.
INDEL_WINDOW = 12


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


def reference_spans(read):
    """
    The read's alignment on the reference as 0-based, half-open intervals:
    aligned bases and deletions, split at each skipped intron (N).
    """
    spans, start, pos = [], read.reference_start, read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 2, 7, 8):  # M, D, =, X
            pos += length
        elif op == 3:  # N
            if pos > start:
                spans.append((start, pos))
            pos += length
            start = pos
    if pos > start:
        spans.append((start, pos))
    return spans


def in_bundle(read, target):
    """Whether a read's reference span, deletions included and introns not,
    overlaps the target's bundle."""
    return read.reference_name == target.chrom and any(
        start < target.hi and end > target.lo for start, end in reference_spans(read))


def read_group(tid, cluster):
    """Name of a cluster's read group (ID and LB) and MI tag: <target_id>_<cluster ID>."""
    return f"{tid}_{cluster}"


def read_group_target(name):
    """The target ID of a read group named by read_group."""
    return name.rsplit("_", 1)[0]


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


##


# PROTOTYPE: position-tolerant indel matching (never merged)
#
# Rule: at an indel target, a read is ALT if it carries an indel of the target's
# type ('+' or '-') whose length is within INDEL_K bp of the target's and whose
# left end lies within `window` bp of the anchor, whatever its sequence. This
# contains genotype-v3's exact rule (offset 0, length within k, prefix agreeing).
# The left end is the anchor base of an insertion (the reference base before it)
# and the first deleted base of a deletion, as bam-readcount places them.
#
# bam-readcount still counts every read at the anchor. For each read that the
# rule makes ALT but v3 does not, tolerant_counts records how bam-readcount saw
# it at the anchor (REF, OTH or not usable), so that count_consensus.py can move
# it to ALT. With window = 0 only anchor entries of the right type and length
# whose sequence disagrees with the alt move (v3: OTH); with window = -1 nothing
# moves.

ANCHOR_REF, ANCHOR_OTH, ANCHOR_NONE, ANCHOR_ALT = "ref", "oth", "none", "alt"
TOLERANT_COLUMNS = ["read_group", "shift_ref", "shift_oth", "shift_none",
                    "v3_ref", "v3_alt", "v3_oth"]


class IndelTarget(NamedTuple):
    """An indel target as the read classifier needs it: p0 is the 0-based
    bam-readcount position (insertion anchor or first deleted base)."""

    target_id: str
    chrom: str
    p0: int
    ref: str
    kind: str  # '+' or '-'
    seq: str   # inserted or deleted bases


def read_indel_targets(bed):
    """Indel targets of a BED, one per target ID, in file order."""
    targets = {}
    with open(bed) as f:
        for line in f:
            if line.strip():
                chrom, start, _, ref, alt = line.rstrip("\n").split("\t")[:5]
                if alt.startswith(("+", "-")):
                    tid = target_id(chrom, start, ref, alt)
                    targets.setdefault(tid, IndelTarget(tid, chrom, int(start) - 1, ref.upper(),
                                                        alt[0], alt[1:].upper()))
    return list(targets.values())


def read_indels(read, p0):
    """
    The read's indels and its base on the 0-based reference position p0.

    Returns:
        Tuple[List[Tuple[str, int, str or None, int]], str]: indels as (kind,
        0-based left end, inserted bases or None for a deletion, length), and
        the read base aligned on p0 (M, = or X), or '' if none is.
    """
    indels, base = [], ""
    seq = read.query_sequence or ""
    rpos, qpos = read.reference_start, 0
    for op, length in read.cigartuples:
        if op in (0, 7, 8):  # M, =, X
            if rpos <= p0 < rpos + length and qpos + p0 - rpos < len(seq):
                base = seq[qpos + p0 - rpos]
            rpos += length
            qpos += length
        elif op == 1:  # I: anchored on the reference base before it
            indels.append(("+", rpos - 1, seq[qpos:qpos + length].upper(), length))
            qpos += length
        elif op == 2:  # D: placed on its first deleted base
            indels.append(("-", rpos, None, length))
            rpos += length
        elif op == 3:  # N
            rpos += length
        elif op == 4:  # S
            qpos += length
    return indels, base


def v3_anchor_class(t, indels, base, k=INDEL_K):
    """
    How bam-readcount and genotype-v3's _indel_allele_counts see one read at the
    target's anchor: ALT, REF, OTH or not usable (NONE).

    An insertion or deletion entry at the anchor decides the class (ALT if it
    matches the alt as in count_consensus.indel_matches, else OTH); an insertion
    entry takes the read off its anchor base. Without one, the anchor base is REF
    if it is the reference base, OTH if another of A, C, G, T, else NONE. A
    deletion at the anchor deletes reference bases, so its sequence is the
    reference and agrees with the alt over their shared part: only its length
    is tested.
    """
    for kind, left, ins_seq, length in indels:
        if left != t.p0:
            continue
        if kind != t.kind or abs(length - len(t.seq)) > k:
            return ANCHOR_OTH
        if kind == "+":
            shared = min(length, len(t.seq))
            return ANCHOR_ALT if ins_seq[:shared] == t.seq[:shared] else ANCHOR_OTH
        return ANCHOR_ALT
    base = base.upper()
    if base == t.ref:
        return ANCHOR_REF
    if base in ("A", "C", "G", "T"):
        return ANCHOR_OTH
    return ANCHOR_NONE


def is_tolerant(t, indels, window=INDEL_WINDOW, k=INDEL_K):
    """Whether the read carries an indel of the target's type, length within k
    of the target's and left end within window bp of the anchor."""
    return any(kind == t.kind and abs(length - len(t.seq)) <= k and abs(left - t.p0) <= window
               for kind, left, _, length in indels)


def tolerant_counts(grouped, indel_targets, window=INDEL_WINDOW, k=INDEL_K):
    """
    Per indel-target cluster: the reads the tolerant rule makes ALT and v3 does
    not, by their anchor class (shift_ref, shift_oth, shift_none), and v3's own
    per-read classes (v3_ref, v3_alt, v3_oth), which should equal bam-readcount's
    r, a and o (a check, not used in the counts).

    Parameters:
        grouped (Iterable[Tuple[str, int, pysam.AlignedSegment]]): group_reads output.
        indel_targets (List[IndelTarget]): see read_indel_targets.

    Returns:
        Dict[str, Dict[str, int]]: read group name -> counts (TOLERANT_COLUMNS[1:]).
    """
    by_id = {t.target_id: t for t in indel_targets}
    out = {}
    for tid, c, read in grouped:
        t = by_id.get(tid)
        if t is None:
            continue
        counts = out.setdefault(read_group(tid, c), dict.fromkeys(TOLERANT_COLUMNS[1:], 0))
        indels, base = read_indels(read, t.p0)
        cls = v3_anchor_class(t, indels, base, k)
        if cls != ANCHOR_NONE:
            counts[f"v3_{cls}"] += 1
        if cls != ANCHOR_ALT and is_tolerant(t, indels, window, k):
            counts[f"shift_{cls}"] += 1
    return out


def write_tolerant_counts(counts, path):
    """Write tolerant_counts as a TSV with TOLERANT_COLUMNS."""
    with open(path, "w") as f:
        f.write("\t".join(TOLERANT_COLUMNS) + "\n")
        for name, c in counts.items():
            f.write("\t".join([name] + [str(c[col]) for col in TOLERANT_COLUMNS[1:]]) + "\n")


##


def write_grouped_bam(bam, bed, output, *, min_mapq, max_ld, ratio, small_reads,
                      tolerant_out=None, indel_window=INDEL_WINDOW):
    """
    Write a cell's grouped reads for bam-readcount -p: one record per read and
    bundle, in a read group per cluster (ID and LB; see read_group), also set
    as the MI tag. The output is coordinate-sorted and indexed. PROTOTYPE: with
    tolerant_out, also write the per-cluster tolerant indel counts there.
    """
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as f:
        header = f.header.to_dict()
        grouped = list(group_reads(f, read_targets(bed), min_mapq=min_mapq, max_ld=max_ld,
                                   ratio=ratio, small_reads=small_reads))
    if tolerant_out is not None:
        write_tolerant_counts(tolerant_counts(grouped, read_indel_targets(bed), indel_window),
                              tolerant_out)

    # Read groups by target, then cluster ID
    order = {tid: i for i, tid in enumerate(dict.fromkeys(tid for tid, _, _ in grouped))}
    names = [read_group(tid, c) for _, c, tid in sorted({(order[tid], c, tid) for tid, c, _ in grouped})]
    cell = Path(bam).name.split(".")[0]
    header["HD"] = {**header.get("HD", {"VN": "1.6"}), "SO": "coordinate"}
    header["RG"] = [{"ID": name, "LB": name, "SM": cell} for name in names]
    out_header = pysam.AlignmentHeader.from_dict(header)

    records = []
    for tid, c, read in grouped:
        record = pysam.AlignedSegment.from_dict(read.to_dict(), out_header)
        record.set_tag("RG", read_group(tid, c), "Z")
        record.set_tag("MI", read_group(tid, c), "Z")
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
    parser.add_argument("--tolerant_out", default=None,
                        help="PROTOTYPE: write per-cluster position-tolerant indel counts here.")
    parser.add_argument("--indel_window", type=int, default=INDEL_WINDOW,
                        help="PROTOTYPE: bp either side of the anchor for a tolerant indel's left end.")
    return parser.parse_args()


def main():
    args = _parse_args()
    write_grouped_bam(args.input, args.bed, args.output, min_mapq=args.min_mapq,
                      max_ld=args.max_ld, ratio=args.ratio, small_reads=args.small_reads,
                      tolerant_out=args.tolerant_out, indel_window=args.indel_window)


if __name__ == "__main__":
    main()
