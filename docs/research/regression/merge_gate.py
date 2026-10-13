"""Merge gate: genotype-v2 fidelity on PALM14E (#41).

v2 and the emulation start from the same tagged on-target PALM14E BAM: the per-cell BAMs that
v2's own SPLIT_BAM wrote in the Nextflow run (--input_type tagged_bam) are read by both.

  G0 tie       v2's own CONSENSUS functions, re-run here per cell (umi_grouping.write_grouped_bam,
               the two bam-readcount passes, count_consensus.count_targets), reproduce the
               pipeline's <sample>.csv WT / MUT / OTH exactly. This makes the per-read cluster
               membership and per-cluster labels below v2's, not a re-implementation.
  a1 grouping  per cell x target, v2's bundle (reads with an MI tag for the target) equals the
               emulation's, and v2's read -> cluster partition equals the emulation's exactly.
               Emulation: primary mapped, MAPQ >= 20, UB-tagged reads in a pysam pileup over the
               target (SNV: the site; indel: target +-k, deletions included, introns not),
               clustered with the measurement helper hybrid.hybrid (LD <= 2, T = 3, downhill).
  a2 caller    the reference spec caller (ladder.call_spec_all) on v2's own WT / MUT / OTH
               reproduces v2's calls and P_s exactly (and p, reported).
  a3 labels    per cluster, v2's class (sub-n_min / WT / MUT / OTH / discarded, from its
               bam-readcount -p blocks) vs the emulation's under rule A on usable reads, computed
               per read: SNV = pileup base, q >= 20, A/C/G/T; indel = the read's site-line entry
               (kit_indel.read_info, no quality filter), ALT by the tolerant match (k = 1). The
               cluster-read-gated helpers (hybrid.consensus, kit_indel.uniform_label) are not used.
  G4 migration v2's calls vs main's on the same subset (palm14e_rerun.csv), for the README note.

Writes AGGREGATE tables (targets as T1..Tn, BED order) to <out>/aggregate and per-cluster /
per-cell tables to <out>/private (HPC only, never committed).
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pysam

import hybrid
from kit_indel import matches, read_info
from ladder import call_spec_all
from measure_umi import close_pairs, levenshtein, load_targets

MAPQ = 20
BQ = 20
NMIN = 3
FMIN = 0.6
K = 1
SMALL = 3
MAXD = 2
LABELLED = ("WT", "MUT", "OTH")
CLASSES = ("sub-n_min", "WT", "MUT", "OTH", "discarded")
V2 = {}  # v2 modules, imported in the workers from --v2-bin


# ---------------------------------------------------------------- v2 side (v2's own code)
def v2_cell(bam, bed, ref, tmp, targets_df):
    """Run v2's CONSENSUS steps on one cell. Returns (tally table, membership, per-cluster r/a/o).

    membership: tid -> {read key -> cluster name}; clusters: (tid, cluster name) -> (r, a, o, class)
    """
    ug, cc = V2["umi_grouping"], V2["count_consensus"]
    grouped = os.path.join(tmp, "grouped.bam")
    ug.write_grouped_bam(bam, bed, grouped, min_mapq=MAPQ, max_ld=MAXD, ratio=2, small_reads=SMALL)
    membership = defaultdict(dict)
    with pysam.AlignmentFile(grouped, "rb") as f:
        for r in f:
            mi = r.get_tag("MI")
            membership[ug.read_group_target(mi)][(r.query_name, r.reference_start)] = mi
    lines = []
    for kind, minbq in (("snv", BQ), ("indel", 0)):
        sub = os.path.join(tmp, f"{kind}.bed")
        with open(bed) as fi, open(sub, "w") as fo:
            for line in fi:
                if line.strip() and (line.split("\t")[4].startswith(("+", "-")) == (kind == "indel")):
                    fo.write(line)
        if os.path.getsize(sub):
            out = subprocess.run(["bam-readcount", "-f", ref, "--site-list", sub, "-b", str(minbq), grouped, "-p"],
                                 check=True, capture_output=True, text=True).stdout
            lines += [l for l in out.splitlines(keepends=True) if l.strip()]
    tally = cc.count_targets(targets_df, lines, FMIN, NMIN)
    clusters = {}
    for line in lines:
        chrom, pos, refb, blocks = cc.parse_readcount_line(line)
        for t in targets_df.itertuples(index=False):
            if (t.chr, t.start, t.ref) != (chrom, pos, refb):
                continue
            tid = V2["genotyping"].target_id(chrom, pos, refb, t.alt)
            for lib, entries in blocks:
                if ug.read_group_target(lib) == tid:
                    r, a, o = cc.allele_counts(entries, refb, t.alt)
                    clusters[(tid, lib)] = (r, a, o, cc.classify_counts(r, a, o, FMIN, NMIN))
    return tally, membership, clusters


# ---------------------------------------------------------------- emulation (independent of v2's code)
def interval(t):
    """0-based, half-open reference interval a read must reach to join the target's bundle."""
    p = t.start - 1
    if t.alt.startswith("-"):
        return p - K, p + len(t.alt) - 1 + K
    if t.alt.startswith("+"):
        return p - K, p + 1 + K
    return p, p + 1


def rule_a(r, a, o):
    u = r + a + o
    if u < NMIN:
        return "sub-n_min"
    for lab, c in (("WT", r), ("MUT", a), ("OTH", o)):
        if c >= NMIN and c / u >= FMIN:
            return lab
    return "discarded"


def indel_allele(entries, t):
    """Allele of one read from its site-line entries (base, then any insertion at the anchor)."""
    if not entries:
        return None
    if entries[0].startswith("-"):
        return "ALT" if matches(entries[0], t.alt, K) else "OTH"
    if len(entries) > 1:  # an insertion read: classified by its insertion, not its anchor base
        return "ALT" if matches(entries[1], t.alt, K) else "OTH"
    b = entries[0]
    return "REF" if b == t.ref else (None if b == "N" else "OTH")


def emu_cell(bam, bai, T, fasta):
    """Emulated bundles, clusters and per-cluster rule-A classes for one cell.

    Returns membership: tid -> {read key -> cluster rep UMI}, clusters: (tid, rep) -> (r, a, o, class)
    """
    kw = dict(truncate=True, stepper="nofilter", ignore_overlaps=False, ignore_orphans=False,
              min_base_quality=0, max_depth=10 ** 9)
    membership, clusters = {}, {}
    with pysam.AlignmentFile(bam, "rb", index_filename=bai) as f:
        for t in T.itertuples(index=False):
            if t.chr not in f.references:
                continue
            lo, hi = interval(t)
            reads, allele = {}, {}
            for col in f.pileup(t.chr, max(0, lo), hi, **kw):
                for pr in col.pileups:
                    a = pr.alignment
                    if (a.flag & 0x904 or a.mapping_quality < MAPQ or not a.has_tag("UB")
                            or pr.is_refskip):
                        continue
                    key = (a.query_name, a.reference_start)
                    reads[key] = a.get_tag("UB")
                    if t.kind == "SNV" and col.reference_pos == t.start - 1 and not pr.is_del \
                            and pr.query_position is not None and a.query_qualities[pr.query_position] >= BQ:
                        b = a.query_sequence[pr.query_position]
                        allele[key] = "REF" if b == t.ref else "ALT" if b == t.alt else "OTH" if b in "ACGT" else None
            if t.kind == "indel":
                w0 = max(0, t.start - 1 - 50)
                seq = fasta.fetch(t.chr, w0, t.start - 1 + 50).upper()
                for r in f.fetch(t.chr, max(0, lo), hi):
                    key = (r.query_name, r.reference_start)
                    if key in reads:
                        allele[key] = indel_allele(read_info(r, t, seq, w0)[0], t)
            if not reads:
                continue
            vc = Counter(reads.values())
            umis = sorted(vc.items(), key=lambda x: (-x[1], x[0]))
            us = [u for u, _ in umis]
            nbrs = defaultdict(list)
            for i, j, _ in sorted(close_pairs(us, MAXD, levenshtein)):
                nbrs[i].append(j)
                nbrs[j].append(i)
            rep, _ = hybrid.hybrid(umis, nbrs, SMALL, True)
            pos = {u: i for i, u in enumerate(us)}
            m = {key: us[rep[pos[ub]]] for key, ub in reads.items()}
            membership[t.tid] = m
            cnt = defaultdict(Counter)
            for key, c in m.items():
                if allele.get(key):
                    cnt[c][allele[key]] += 1
            for c in set(m.values()):
                r, a, o = cnt[c]["REF"], cnt[c]["ALT"], cnt[c]["OTH"]
                clusters[(t.tid, c)] = (r, a, o, rule_a(r, a, o))
    return membership, clusters


# ---------------------------------------------------------------- per cell comparison
def partition(m):
    p = defaultdict(set)
    for key, c in m.items():
        p[c].add(key)
    return {frozenset(v): c for v, c in p.items()}


def compare_cell(args):
    cell, bam, bed, ref, T, tmp_root = args
    targets_df = T[["chr", "start", "end", "ref", "alt", "gene"]]
    tmp = tempfile.mkdtemp(dir=tmp_root)
    try:
        tally, v2m, v2c = v2_cell(bam, bed, ref, tmp, targets_df)
        bai = os.path.join(tmp, "cell.bai")  # SPLIT_BAM writes no index; the work dir is left untouched
        pysam.index(bam, bai)
        emm, emc = emu_cell(bam, bai, T, pysam.FastaFile(ref))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    tally = tally.assign(cell=cell)
    tid_v2 = dict(zip(T.tid, [V2["genotyping"].target_id(c, s, r, a) for c, s, r, a in zip(T.chr, T.start, T.ref, T.alt)]))
    stats, private = [], []
    for t in T.itertuples(index=False):
        vm, em = v2m.get(tid_v2[t.tid], {}), emm.get(t.tid, {})
        vp, ep = partition(vm), partition(em)
        same_bundle = set(vm) == set(em)
        same_part = set(vp) == set(ep)
        s = {"cell": cell, "tid": t.tid, "reads_v2": len(vm), "reads_emu": len(em),
             "bundle_equal": same_bundle, "partition_equal": same_part,
             "clusters_v2": len(vp), "clusters_emu": len(ep),
             "clusters_matched": len(set(vp) & set(ep))}
        pairs = Counter()
        for members in set(vp) | set(ep):
            vc_, ec_ = vp.get(members), ep.get(members)
            v = v2c.get((tid_v2[t.tid], vc_), (0, 0, 0, "sub-n_min")) if vc_ is not None else None
            e = emc[(t.tid, ec_)] if ec_ is not None else None
            if v is None or e is None:
                pairs[("unmatched_v2" if e is None else "unmatched_emu", "")] += 1
                private.append({"cell": cell, "tid": t.tid, "reads": len(members), "v2": v, "emu": e, "why": "partition"})
                continue
            pairs[(v[3], e[3])] += 1
            pairs[("rao_equal", "")] += v[:3] == e[:3]
            if v[3] != e[3] or v[:3] != e[:3]:
                private.append({"cell": cell, "tid": t.tid, "reads": len(members), "v2": v, "emu": e,
                                "why": "label" if v[3] != e[3] else "counts"})
        s["pairs"] = dict(pairs)
        stats.append(s)
    return tally, stats, private


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell-bams", required=True, help="directory with v2's SPLIT_BAM per-cell BAMs")
    ap.add_argument("--v2-bin", required=True, help="genotype-v2 bin/long_reads")
    ap.add_argument("--v2-results", required=True, help="v2 <outdir>/<sample> with <sample>.csv and _targets.csv")
    ap.add_argument("--sample", default="palm14e")
    ap.add_argument("--bed", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--main-calls", required=True, help="main's calls on the same subset (palm14e_rerun.csv)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--max-cells", type=int, default=0)
    a = ap.parse_args()
    sys.path.insert(0, a.v2_bin)
    agg, prv = os.path.join(a.out, "aggregate"), os.path.join(a.out, "private")
    os.makedirs(agg, exist_ok=True)
    os.makedirs(prv, exist_ok=True)
    T = load_targets(a.bed)
    bams = sorted(glob.glob(os.path.join(a.cell_bams, "*.bam")))
    if a.max_cells:
        bams = bams[:a.max_cells]
    cells = [os.path.basename(b).split(".")[0] for b in bams]
    tmp_root = tempfile.mkdtemp(prefix="gate_", dir=os.environ.get("TMPDIR", "/tmp"))
    print("cells", len(cells), flush=True)

    jobs = [(c, b, a.bed, a.ref, T, tmp_root) for c, b in zip(cells, bams)]
    with Pool(a.threads, initializer=_init, initargs=(a.v2_bin,)) as pool:
        out = []
        for i, res in enumerate(pool.imap(compare_cell, jobs, chunksize=4)):
            out.append(res)
            if (i + 1) % 500 == 0:
                print("cells done", i + 1, flush=True)
    shutil.rmtree(tmp_root, ignore_errors=True)
    _init(a.v2_bin)
    tally = pd.concat([o[0] for o in out], ignore_index=True)
    stats = pd.DataFrame([s for o in out for s in o[1]])
    private = pd.DataFrame([p for o in out for p in o[2]])
    tid_of = {V2["genotyping"].target_id(c, s, r, x): t for c, s, r, x, t in zip(T.chr, T.start, T.ref, T.alt, T.tid)}
    kind_of = dict(zip(T.tid, T.kind))

    # ---------------- G0 tie: harness tallies vs the pipeline's <sample>.csv
    v2 = pd.read_csv(os.path.join(a.v2_results, f"{a.sample}.csv"))
    v2["tid"] = v2.target_id.map(tid_of)
    v2 = v2[v2.cell.isin(set(cells))]
    tally["tid"] = [tid_of[V2["genotyping"].target_id(c, s, r, x)] for c, s, r, x in
                    zip(tally.chr, tally.start, tally.ref, tally.alt)]
    j = v2.set_index(["cell", "tid"])[list(LABELLED)].join(
        tally.set_index(["cell", "tid"])[list(LABELLED)], rsuffix="_h", how="outer")
    g0 = []
    for tid in list(T.tid) + ["all"]:
        s = j if tid == "all" else j[j.index.get_level_values("tid") == tid]
        eq = (s[list(LABELLED)].values == s[[c + "_h" for c in LABELLED]].values).all(axis=1)
        g0.append({"tid": tid, "rows": len(s), "rows_identical": int(eq.sum()), "rows_differing": int((~eq).sum()),
                   **{f"{c}_pipeline": int(s[c].sum()) for c in LABELLED},
                   **{f"{c}_harness": int(s[c + "_h"].sum()) for c in LABELLED}})
    pd.DataFrame(g0).to_csv(os.path.join(agg, "G0_tie.csv"), index=False)

    # ---------------- a1 grouping, a3 labels per target
    a1, a3, conf = [], [], []
    for tid in T.tid:
        s = stats[stats.tid == tid]
        cov = s[(s.reads_v2 > 0) | (s.reads_emu > 0)]
        a1.append({"tid": tid, "kind": kind_of[tid], "cell_bundles": len(cov),
                   "reads_v2": int(s.reads_v2.sum()), "reads_emu": int(s.reads_emu.sum()),
                   "bundles_equal": int(cov.bundle_equal.sum()), "partitions_equal": int(cov.partition_equal.sum()),
                   "clusters_v2": int(s.clusters_v2.sum()), "clusters_emu": int(s.clusters_emu.sum()),
                   "clusters_matched": int(s.clusters_matched.sum()),
                   "a1_pass": bool(cov.partition_equal.all())})
        P = Counter()
        for p in s.pairs:
            P.update(p)
        matched = sum(n for (v, e), n in P.items() if v in CLASSES)
        agree = sum(n for (v, e), n in P.items() if v in CLASSES and v == e)
        lab = sum(n for (v, e), n in P.items() if v in CLASSES and (v in LABELLED or e in LABELLED))
        lab_agree = sum(n for (v, e), n in P.items() if v in CLASSES and v == e and v in LABELLED)
        a3.append({"tid": tid, "kind": kind_of[tid], "clusters_matched": matched,
                   "class_agree": agree, "class_agree_frac": agree / matched if matched else np.nan,
                   "clusters_labelled_either": lab, "label_agree": lab_agree,
                   "label_agree_frac": lab_agree / lab if lab else np.nan,
                   "rao_equal": int(P[("rao_equal", "")]),
                   "a3_pass": (lab_agree / lab if lab else 1.0) >= 0.99})
        for (v, e), n in P.items():
            if v in CLASSES:
                conf.append({"tid": tid, "v2": v, "emu": e, "clusters": int(n)})
    pd.DataFrame(a1).to_csv(os.path.join(agg, "A1_grouping.csv"), index=False)
    pd.DataFrame(a3).to_csv(os.path.join(agg, "A3_labels.csv"), index=False)
    pd.DataFrame(conf).to_csv(os.path.join(agg, "A3_confusion.csv"), index=False)
    if len(private):
        private.to_csv(os.path.join(prv, "cluster_disagreements.csv.gz"), index=False)
    stats.drop(columns="pairs").to_csv(os.path.join(prv, "cell_target_grouping.csv.gz"), index=False)

    # ---------------- a2 caller: reference spec caller on v2's own counts
    tab = v2[["cell", "tid"] + list(LABELLED)].copy()
    ref_calls, ps = call_spec_all(tab, T)
    r2 = v2.set_index(["cell", "tid"]).join(ref_calls.set_index(["cell", "tid"])[["p", "genotype"]], rsuffix="_ref")
    gt_v2 = r2.genotype.fillna("NA")
    ps_ref = r2.index.get_level_values("tid").map(lambda t: ps[t][1])
    cov = (r2.MUT + r2.WT) > 0
    targets = pd.read_csv(os.path.join(a.v2_results, f"{a.sample}_targets.csv"))
    targets["tid"] = targets.target_id.map(tid_of)
    a2 = []
    for tid in list(T.tid) + ["all"]:
        m = np.ones(len(r2), bool) if tid == "all" else (r2.index.get_level_values("tid") == tid)
        row = {"tid": tid, "rows": int(m.sum()),
               "calls_equal": int((gt_v2[m] == r2.genotype_ref[m]).sum()),
               "P_s_equal": int((r2.P_s.values[m] == np.asarray(ps_ref)[m]).sum()),
               "p_equal_covered": int((r2.p[m & cov] == r2.p_ref[m & cov]).sum()),
               "rows_covered": int((m & cov).sum()),
               "p_max_abs_diff": float((r2.p[m & cov] - r2.p_ref[m & cov]).abs().max()) if (m & cov).any() else 0.0}
        if tid != "all":
            tr = targets[targets.tid == tid].iloc[0]
            row.update({"N": int(tr.N), "w": int(tr.w), "P_s_v2_targets": tr.P_s, "P_s_ref": ps[tid][1],
                        "q_hat_ref": ps[tid][0]})
        row["a2_pass"] = row["calls_equal"] == row["rows"] and row["P_s_equal"] == row["rows"]
        a2.append(row)
    pd.DataFrame(a2).to_csv(os.path.join(agg, "A2_caller.csv"), index=False)

    # ---------------- G4 migration numbers: v2 vs main on the same subset
    mc = pd.read_csv(a.main_calls).rename(columns={"barcode": "cell"})
    tmap = {(r.chr, r.start, r.alt): r.tid for r in T.itertuples()}
    mc["tid"] = [tmap.get((c, s, x)) for c, s, x in zip(mc.chr, mc.start, mc.alt)]
    mc = mc.dropna(subset=["tid"])
    mj = mc.set_index(["cell", "tid"]).genotype.fillna("NA").rename("main").to_frame().join(
        gt_v2.rename("v2"), how="inner")
    changed = mj.main != mj.v2
    hard = mj.main.isin(["MUT", "WT"]) & mj.v2.isin(["MUT", "WT"]) & changed
    pc = changed.groupby(level="cell").any()
    g4 = [{"metric": "cells (common grid)", "value": int(pc.size)},
          {"metric": "cells in v2 / main tables", "value": f"{v2.cell.nunique()} / {mc.cell.nunique()}"},
          {"metric": "MUT calls main", "value": int((mj.main == "MUT").sum())},
          {"metric": "MUT calls v2", "value": int((mj.v2 == "MUT").sum())},
          {"metric": "WT calls main", "value": int((mj.main == "WT").sum())},
          {"metric": "WT calls v2", "value": int((mj.v2 == "WT").sum())},
          {"metric": "cells changing >= 1 call", "value": int(pc.sum())},
          {"metric": "hard flips MUT>WT", "value": int((hard & (mj.main == "MUT")).sum())},
          {"metric": "hard flips WT>MUT", "value": int((hard & (mj.main == "WT")).sum())},
          {"metric": "v2 MUT calls with m = 1", "value": int(((v2.genotype == "MUT") & (v2.MUT == 1)).sum())}]
    for tid in T.tid:
        s = mj[mj.index.get_level_values("tid") == tid]
        g4.append({"metric": f"{tid} MUT calls main -> v2", "value": f"{int((s.main == 'MUT').sum())} -> {int((s.v2 == 'MUT').sum())}"})
    pd.DataFrame(g4).to_csv(os.path.join(agg, "G4_migration.csv"), index=False)

    # ---------------- verdict
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    for name in ["G0_tie", "A1_grouping", "A2_caller", "A3_labels", "G4_migration"]:
        print(f"--- {name}")
        print(pd.read_csv(os.path.join(agg, f"{name}.csv")).to_string(index=False))
    print("done", flush=True)


def _init(v2_bin):
    if v2_bin not in sys.path:
        sys.path.insert(0, v2_bin)
    import count_consensus
    import genotyping
    import umi_grouping
    V2.update(umi_grouping=umi_grouping, count_consensus=count_consensus, genotyping=genotyping)


if __name__ == "__main__":
    main()
