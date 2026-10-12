#!/usr/bin/env python
"""UMI grouping within one bundle (cell x target): the Levenshtein hybrid merge.

UMI a absorbs b if LD(a, b) <= max_ld and either
  n_a >= ratio * n_b - 1                (directional ratio), or
  n_b <= small_reads and n_a >= n_b     (small clause, kept downhill)
where n is the number of reads per raw UB. Clusters grow from the most abundant
UMI through its absorbed UMIs; ties break by count, then by UMI sequence. The
allele is never used. Measured on PALM14E in docs/research/umi-hybrid-rule.md.
"""

from collections import defaultdict
from itertools import combinations


def levenshtein(a, b, max_d):
    """Levenshtein distance between a and b, or max_d + 1 if it exceeds max_d.

    Banded DP: only cells within max_d of the diagonal are filled, and the scan
    stops as soon as a whole row exceeds max_d.
    """
    if a == b:
        return 0
    if abs(len(a) - len(b)) > max_d:
        return max_d + 1
    out = max_d + 1
    # prev[j] is the distance between a[:i-1] and b[:j]; cells outside the band are out
    prev = [j if j <= max_d else out for j in range(len(b) + 1)]
    for i in range(1, len(a) + 1):
        lo, hi = max(1, i - max_d), min(len(b), i + max_d)
        cur = [out] * (len(b) + 1)
        if i <= max_d:
            cur[0] = i
        for j in range(lo, hi + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1]), out)
        if min(cur[max(0, lo - 1) : hi + 1]) > max_d:
            return out
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
        raise ValueError("UMI read counts must be positive")
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
