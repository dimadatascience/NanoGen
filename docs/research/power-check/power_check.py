#!/usr/bin/env python3
"""Outcome-blind power check (map #72, ticket #76).

Coverage only: no allele is ever read. Per cell x target, n_proxy is the number of
genotype-v3 UMI clusters (umi_grouping.py imported from genotype-v3@61d8809) with
>= min_reads usable reads at the site. A read is usable at an SNV if it has an
aligned base with quality >= snv_min_base_quality or a deletion over the site;
at an indel every bundle read is usable (v3 runs bam-readcount -b 0 there).
This is an upper bound on v3's n = WT + MUT: it also counts the clusters that the
consensus step labels OTH or discards (no allele reaches min_fraction). The gap
is measured on PALM22E/PALM09E, whose v3 runs are already scored.

  coverage   BAM -> per-cell table (barcodes; stays on the HPC)
  summarise  per-cell table -> aggregate per-target table (no barcodes), with the
             pre-registered evaluability (#64) and calibration power under a fixed
             grid of contamination scenarios; with --v3-csv also the proxy-vs-v3 gap.

GENOTYPING_DIR must point at genotype-v3@61d8809 bin/long_reads (imports umi_grouping,
genotyping.call via score_lib).
"""
import argparse
import csv
import gzip
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.environ["GENOTYPING_DIR"])

BARE = re.compile(r"[ACGTN]{16}")
# Contamination scenarios (θ = P_s + c): fixed before any new-sample look, spanning the
# c_s seen on the last map's targets (0.09 %-4 %). P_s at the v3 floor.
C_GRID = [0.0025, 0.005, 0.01, 0.02, 0.04]
P_S = 1e-4
ALPHA_CAL, RATIO_CAP = 0.05, 1.5  # criterion 1 of #64: pass iff p_cal >= 0.05 and O/E <= 1.5
POWER = 0.8


def bare(c):
    m = BARE.search(c)
    return m.group(0) if m else None


def read_list(path):
    """Barcodes of a one-column list or of a CSV with a 'cell' / 'cell16' column (or no header)."""
    out = set()
    with open(path) as f:
        rows = list(csv.reader(f))
    if not rows:
        return out
    head = [h.strip().strip('"') for h in rows[0]]
    col = next((head.index(k) for k in ("cell", "cell16") if k in head), 0)
    for r in rows[1:] if any(k in head for k in ("cell", "cell16")) else rows:
        if len(r) > col and bare(r[col]):
            out.add(bare(r[col]))
    return out


def read_labels(path):
    """barcode -> ct from a CSV with cell16,ct columns, or a header-less cell,ct file."""
    lab = {}
    with open(path) as f:
        rows = list(csv.reader(f))
    head = [h.strip().strip('"') for h in rows[0]]
    if "ct" in head:
        ic, it = head.index("cell16"), head.index("ct")
        body = rows[1:]
    else:
        ic, it, body = 0, 1, rows
    for r in body:
        b = bare(r[ic])
        if b:
            lab[b] = r[it]
    return lab


# ---------------------------------------------------------------- coverage

def site_usable(read, pos, min_bq):
    """Whether the read has a usable observation at reference pos: an aligned base with
    quality >= min_bq, or a deletion over pos. The base itself is never inspected."""
    ref, q = read.reference_start, 0
    quals = read.query_qualities
    for op, ln in read.cigartuples:
        if op in (0, 7, 8):
            if ref <= pos < ref + ln:
                return quals is None or quals[q + pos - ref] >= min_bq
            ref += ln
            q += ln
        elif op in (1, 4):
            q += ln
        elif op == 2:
            if ref <= pos < ref + ln:
                return True
            ref += ln
        elif op == 3:
            if ref <= pos < ref + ln:
                return False
            ref += ln
    return False


def coverage(a):
    import pysam
    from umi_grouping import cluster_umis, in_bundle, is_primary_mapped, read_targets

    universe = read_list(a.qc) | read_list(a.negatives)
    targets = read_targets(a.bed)
    bam = pysam.AlignmentFile(a.bam)
    opener = gzip.open if a.out.endswith(".gz") else open
    with opener(a.out, "wt") as out:
        out.write("T\ttarget_id\tcell\tn_proxy\tclusters\tbundle_reads\n")
        for i, t in enumerate(targets):
            snv = not t.target_id.rsplit(":", 1)[1].split(">")[1][:1] in "+-"
            members = defaultdict(list)  # cell -> [(UB, usable)]
            for r in bam.fetch(t.chrom, max(0, t.lo), t.hi):
                if not is_primary_mapped(r, a.min_mapq) or not r.has_tag("UB") or not r.has_tag("CB"):
                    continue
                cb = bare(r.get_tag("CB"))
                if cb not in universe or not in_bundle(r, t):
                    continue
                use = site_usable(r, t.lo, a.snv_min_bq) if snv else True
                members[cb].append((r.get_tag("UB"), use))
            for cb, m in members.items():
                cl = cluster_umis(Counter(u for u, _ in m), max_ld=a.max_ld, ratio=a.ratio,
                                  small_reads=a.small_reads)
                usable = Counter(cl[u] for u, ok in m if ok)
                n = sum(1 for v in usable.values() if v >= a.min_reads)
                out.write(f"T{i + 1}\t{t.target_id}\t{cb}\t{n}\t{len(set(cl.values()))}\t{len(m)}\n")
            print(f"T{i + 1} done: {len(members)} cells", flush=True)


# ---------------------------------------------------------------- power

def crit_k(E):
    """Smallest O with Poisson tail P(X >= O | E) < ALPHA_CAL."""
    from scipy.stats import poisson
    k = int(np.floor(E))
    while poisson.sf(k - 1, E) >= ALPHA_CAL:
        k += 1
    return k


def power_numbers(E):
    """Calibration power at null expectation E (O ~ Poisson(r E) as a stand-in for the
    Poisson-binomial tail of #64): r80 (smallest O/E detected by p_cal < 0.05 with 80 %
    power), the criterion-1 fail rate (p_cal < 0.05 or O/E > 1.5) at r = 1 and r = 1.5,
    and r80_lo (largest O/E < 1 detected by the lower tail with 80 % power)."""
    from scipy.stats import poisson
    if E <= 0:
        return dict(r80=np.inf, fail_r1=np.nan, fail_r15=np.nan, r80_lo=0.0)
    k = crit_k(E)
    k_fail = min(k, int(np.floor(RATIO_CAP * E)) + 1)
    r80 = next((r for r in np.arange(1.0, 20.0, 0.01) if poisson.sf(k - 1, r * E) >= POWER), np.inf)
    lo = int(poisson.ppf(ALPHA_CAL, E))  # P(X <= lo-1) < alpha
    while lo >= 0 and poisson.cdf(lo, E) >= ALPHA_CAL:
        lo -= 1
    r80_lo = (next((r for r in np.arange(0.99, 0.0, -0.01) if poisson.cdf(lo, r * E) >= POWER), 0.0)
              if lo >= 0 else 0.0)
    return dict(r80=float(r80), fail_r1=float(poisson.sf(k_fail - 1, E)),
                fail_r15=float(poisson.sf(k_fail - 1, RATIO_CAP * E)), r80_lo=float(r80_lo))


def e_star():
    """Smallest null expectation E at which a 1.5x miscalibration is detected (p_cal < 0.05)
    with 80 % power."""
    from scipy.stats import poisson
    for E in np.arange(0.5, 200, 0.1):
        if poisson.sf(crit_k(E) - 1, RATIO_CAP * E) >= POWER:
            return float(E)
    return np.inf


def c_min(Es, es):
    """Smallest c with E_B(c) >= es, log-log interpolated over C_GRID; inf above the grid."""
    for j, e in enumerate(Es):
        if e >= es:
            if j == 0 or Es[j - 1] <= 0:
                return C_GRID[j]
            return float(np.exp(np.interp(np.log(es), np.log([Es[j - 1], e]),
                                          np.log([C_GRID[j - 1], C_GRID[j]]))))
    return np.inf


def summarise(a):
    import score_lib as L  # imports genotyping.call (genotype-v3@61d8809)

    qc, neg_all = read_list(a.qc), read_list(a.negatives)
    neg = sorted(neg_all & qc)
    lab = read_labels(a.labels)
    leuk = {b for b, c in lab.items() if c == "Leukemic"} & qc
    splits = [L.split(neg, s) for s in L.SEEDS]
    cov = defaultdict(dict)
    order = {}
    with (gzip.open if a.cells.endswith(".gz") else open)(a.cells, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["cell"] in qc:
                cov[r["T"]][r["cell"]] = int(r["n_proxy"])
            order.setdefault(r["T"], r["target_id"])
    v3 = None
    if a.v3_csv:
        v3 = defaultdict(dict)
        with open(a.v3_csv) as f:
            for r in csv.DictReader(f):
                b = bare(r["cell"])
                if b in qc:
                    v3[r["target_id"]][b] = int(r["WT"]) + int(r["MUT"])  # MUT never reported alone
    es = e_star()
    rows = []
    targets = [f"T{i}" for i in range(1, a.n_targets + 1)]
    for T in targets:
        n = cov.get(T, {})
        tid = order.get(T, "")
        row = dict(sample=a.sample, target=T, neg_qc=len(neg), leuk_qc=len(leuk))

        def evaluate(nm, scale_u=1.0, scale_c=1.0, tag=""):
            ev, au, bc, eb = [], [], [], {c: [] for c in C_GRID}
            for A, B in splits:
                nA = sum(nm.get(b, 0) for b in A) * scale_u
                nb = [nm[b] for b in B if nm.get(b, 0) >= 1]
                ev.append(nA >= L.MIN_A_UMIS and len(nb) * scale_c >= L.MIN_B_CELLS)
                au.append(nA)
                bc.append(len(nb) * scale_c)
                for c in C_GRID:
                    th = min(1.0, P_S + c)
                    eb[c].append(float(L.p_mut(np.array(nb, int), th, 1, th).sum()) * scale_u if nb else 0.0)
            d = {f"eval_seeds{tag}": int(sum(ev)), f"evaluable{tag}": int(sum(ev)) >= L.MIN_EVAL_SEEDS,
                 f"A_umis_med{tag}": float(np.median(au)), f"B_cells_med{tag}": float(np.median(bc))}
            for c in C_GRID:
                d[f"E_B_c{c}{tag}"] = float(np.median(eb[c]))
            return d

        row.update(evaluate(n))
        neg_cov = [n[b] for b in neg if n.get(b, 0) >= 1]
        row.update(neg_covered=len(neg_cov), neg_umis=int(sum(neg_cov)),
                   leuk_covered=sum(1 for b in leuk if n.get(b, 0) >= 1),
                   leuk_umis=int(sum(n.get(b, 0) for b in leuk)))
        if a.scale_umis != 1.0 or a.scale_cells != 1.0:
            row.update(evaluate(n, a.scale_umis, a.scale_cells, "_cal"))
            row["leuk_covered_cal"] = row["leuk_covered"] * a.scale_cells
        tag = "_cal" if f"E_B_c{C_GRID[0]}_cal" in row else ""
        # smallest c at which E_B reaches E* (1.5x detectable with 80 % power)
        row["c_min_1.5x"] = c_min([row[f"E_B_c{c}{tag}"] for c in C_GRID], es)
        for c in C_GRID:
            for k, v in power_numbers(row[f"E_B_c{c}{tag}"]).items():
                row[f"{k}_c{c}"] = v
        if v3 is not None:
            m = v3.get(tid, {})
            nv = [m[b] for b in neg if m.get(b, 0) >= 1]
            row.update(v3_neg_covered=len(nv), v3_neg_umis=int(sum(nv)),
                       v3_leuk_covered=sum(1 for b in leuk if m.get(b, 0) >= 1),
                       v3_leuk_umis=int(sum(m.get(b, 0) for b in leuk)))
            both = [b for b in qc if n.get(b, 0) >= 1 or m.get(b, 0) >= 1]
            row["cells_exact"] = float(np.mean([n.get(b, 0) == m.get(b, 0) for b in both])) if both else np.nan
            row["cells_v3_gt_proxy"] = int(sum(m.get(b, 0) > n.get(b, 0) for b in both))
            ev3 = evaluate(m, tag="_v3")
            row.update({k: ev3[k] for k in ("eval_seeds_v3", "evaluable_v3", "A_umis_med_v3", "B_cells_med_v3")})
        rows.append(row)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(a.out, "w") as f:
        w = csv.DictWriter(f, keys, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"E* = {es:.1f} expected null calls in B for 80 % power at O/E 1.5")


def main():
    p = argparse.ArgumentParser()
    sp = p.add_subparsers(dest="cmd", required=True)
    c = sp.add_parser("coverage")
    for k in ("bam", "bed", "qc", "negatives", "out"):
        c.add_argument(f"--{k}", required=True)
    c.add_argument("--min-mapq", type=int, default=20)
    c.add_argument("--max-ld", type=int, default=2)
    c.add_argument("--ratio", type=float, default=2)
    c.add_argument("--small-reads", type=int, default=3)
    c.add_argument("--min-reads", type=int, default=3)
    c.add_argument("--snv-min-bq", type=int, default=20)
    s = sp.add_parser("summarise")
    for k in ("sample", "cells", "qc", "negatives", "labels", "out"):
        s.add_argument(f"--{k}", required=True)
    s.add_argument("--n-targets", type=int, required=True)
    s.add_argument("--v3-csv")
    s.add_argument("--scale-umis", type=float, default=1.0)
    s.add_argument("--scale-cells", type=float, default=1.0)
    a = p.parse_args()
    a.min_mapq = getattr(a, "min_mapq", 20)
    {"coverage": coverage, "summarise": summarise}[a.cmd](a)


if __name__ == "__main__":
    main()
