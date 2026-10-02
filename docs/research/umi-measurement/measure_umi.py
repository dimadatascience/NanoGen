"""Measure UMI splitting and window double-counting on PALM14E (wayfinder ticket #9).

Reads the per-cell intermediates kept by stage1_consensus.sh (cells/<CB>/filtered_input.bam,
grouped_reads.tsv, allcounts.count) and writes AGGREGATE tables only. Targets are
reported as T1..Tn (BED order); the mapping to coordinates stays on the HPC.

Checks (research doc §4):
  1 umi_tools bundle sizes, BX != UB rate
  2 splitting of (cell, target, exact UB) across windows / UG; read start/end spread
  3 UMI distance spectrum within (cell, target) vs cross-cell null
  4 allele concordance of candidate merges
  5 what n_min removes, UMI counts and calls under current / B / B+C grouping
  6 secondary + supplementary records
plus V: emulation check of per-group allele counts against bam-readcount.
"""
import argparse
import glob
import os
import random
import re
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np
import pandas as pd
import pysam
from scipy.stats import nbinom

MAPQ = 20
BQ = 20
WIN = 100
FMIN = 0.6
ERR_FLOOR = 1e-4
NMINS = (1, 2, 3)


# ---------------------------------------------------------------- utils
def levenshtein(a, b, cap=4):
    if a == b:
        return 0
    if abs(len(a) - len(b)) >= cap:
        return cap
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        for j, cb in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
        prev = cur
    return min(prev[-1], cap)


def hamming(a, b):
    if len(a) != len(b):
        return 99
    return sum(x != y for x, y in zip(a, b))


def load_targets(bed):
    t = pd.read_csv(bed, sep="\t", header=None, names=["chr", "start", "end", "ref", "alt", "gene"])
    t["tid"] = [f"T{i+1}" for i in range(len(t))]
    t["kind"] = np.where(t.alt.str.startswith(("-", "+")), "indel", "SNV")
    return t


def site_alleles(f, t):
    """Yield (alignment, allele) for every record spanning target t, emulating bam-readcount
    at the site via pileup. allele in REF/ALT/OTH or None (not counted: low BQ, N, ref-skip,
    deleted base)."""
    p0 = t.start - 1
    kw = dict(truncate=True, stepper="nofilter", ignore_overlaps=False, ignore_orphans=False,
              min_base_quality=0, max_depth=100000000)
    ev = defaultdict(list)
    if t.kind == "indel":
        for col in f.pileup(t.chr, max(0, p0 - 1), p0 + 2, **kw):
            for pr in col.pileups:
                if pr.indel != 0:
                    a = pr.alignment
                    ev[(a.query_name, a.flag, a.reference_start)].append(pr.indel)
    alen = len(t.alt) - 1
    sign = -1 if t.alt[0] == "-" else 1
    for col in f.pileup(t.chr, p0, p0 + 1, **kw):
        for pr in col.pileups:
            a = pr.alignment
            al = None
            if t.kind == "indel":
                e = ev.get((a.query_name, a.flag, a.reference_start), [])
                if any(x * sign > 0 and abs(abs(x) - alen) <= 1 for x in e):
                    al = "ALT"
                elif e:
                    al = "OTH"
            if al is None and not pr.is_del and not pr.is_refskip and pr.query_position is not None:
                qp = pr.query_position
                if a.query_qualities[qp] >= BQ:
                    b = a.query_sequence[qp]
                    if b != "N":
                        if b == t.ref:
                            al = "REF"
                        elif t.kind == "SNV" and b == t.alt:
                            al = "ALT"
                        else:
                            al = "OTH"
            yield a, al


def close_pairs(us, k, dist):
    """Pairs (i, j, d) of distinct strings in us with dist <= k, via deletion neighbourhoods
    (any pair within Levenshtein/Hamming k shares a variant with <= k deletions each)."""
    idx = defaultdict(list)
    for i, u in enumerate(us):
        vs = {u}
        for n in range(1, k + 1):
            for pos in combinations(range(len(u)), n):
                vs.add("".join(c for m, c in enumerate(u) if m not in pos))
        for v in vs:
            idx[v].append(i)
    cand = set()
    for lst in idx.values():
        if len(lst) > 1:
            for x in range(len(lst)):
                for y in range(x + 1, len(lst)):
                    cand.add((lst[x], lst[y]))
    out = []
    for i, j in cand:
        d = dist(us[i], us[j])
        if d <= k:
            out.append((i, j, d))
    return out


def label(counts, nmin, fmin=FMIN):
    """Uniform consensus rule (spec): n >= nmin and top allele fraction >= fmin."""
    n = sum(counts.values())
    if n < nmin or n == 0:
        return None
    a, c = max(counts.items(), key=lambda x: x[1])
    if c / n < fmin:
        return None
    return {"REF": "WT", "ALT": "MUT", "OTH": "OTH"}[a]


def directional(umis, maxd, dist):
    """umi_tools-style directional clustering: a -> b if d(a,b)<=maxd and n_a >= 2 n_b - 1.
    umis: dict umi -> count. Returns dict umi -> representative."""
    us = sorted(umis, key=lambda u: (-umis[u], u))
    nb = defaultdict(list)
    for i, j, _ in close_pairs(us, maxd, dist):
        nb[i].append(j)
        nb[j].append(i)
    rep = {}
    for i, u in enumerate(us):
        if u in rep:
            continue
        rep[u] = u
        stack = [i]
        while stack:
            x = stack.pop()
            for y in nb[x]:
                v = us[y]
                if v in rep:
                    continue
                if umis[us[x]] >= 2 * umis[v] - 1:
                    rep[v] = u
                    stack.append(y)
    return rep


def call_main(df):
    """main@c442a32 genotyping.py rule on a WT/MUT/OTH table (per-cell err)."""
    g = df.groupby("cell")
    err = (g["OTH"].sum() / (g[["MUT", "WT", "OTH"]].sum().sum(axis=1) + 1e-8) / 2).rename("err")
    df = df.merge(err, left_on="cell", right_index=True)
    mn = df.loc[df.err > 0, "err"].min() if (df.err > 0).any() else ERR_FLOOR
    df.loc[df.err == 0, "err"] = mn
    e = np.maximum(df.err.values, ERR_FLOOR)
    wt, mut = df.WT.values, df.MUT.values
    with np.errstate(all="ignore"):
        p = nbinom.cdf(k=wt, n=mut, p=e)
    gt = np.full(p.size, "NA", dtype=object)
    gt[p < 0.1] = "MUT"
    gt[np.isnan(p) & (wt >= 2)] = "WT"
    gt[p >= 0.9] = "WT"
    gt[(wt == 0) & (mut == 0)] = "NA"
    df["genotype"] = gt
    return df


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--null-pairs", type=int, default=200000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    random.seed(a.seed)
    os.makedirs(a.out, exist_ok=True)
    T = load_targets(a.bed)
    targets = list(T.itertuples(index=False))

    rows = []  # one per (record, target covered)
    bundles = Counter()
    bundle_umis = defaultdict(set)
    n_grouped = n_bx_ne_ub = 0
    rc_groups = {}  # (cell, lb, tid) -> bam-readcount counts REF/ALT/OTH
    rec_stats = Counter()

    celldirs = sorted(glob.glob(os.path.join(a.cells, "*")))
    for cd in celldirs:
        cell = os.path.basename(cd)
        bam = os.path.join(cd, "filtered_input.bam")
        if not os.path.exists(bam):
            continue
        # check 1: umi_tools group-out
        gpath = os.path.join(cd, "grouped_reads.tsv")
        if os.path.exists(gpath):
            g = pd.read_csv(gpath, sep="\t")
            for (ctg, pos), sub in g.groupby(["contig", "position"]):
                bundles[(cell, ctg, pos)] = len(sub)
                bundle_umis[(cell, ctg, pos)] = set(sub.umi)
            n_grouped += len(g)
            n_bx_ne_ub += int((g.umi != g.final_umi).sum())
        # per-record table (needs a sorted, indexed BAM)
        if not os.path.exists(bam + ".bai"):
            pysam.sort("-o", os.path.join(cd, "sorted_fi.bam"), bam)
            bam = os.path.join(cd, "sorted_fi.bam")
            pysam.index(bam)
        with pysam.AlignmentFile(bam, "rb") as f:
            for r in f.fetch(until_eof=True):
                if r.is_unmapped:
                    continue
                kind = "supp" if r.is_supplementary else ("sec" if r.is_secondary else "prim")
                rec_stats[(kind, r.mapping_quality >= MAPQ)] += 1
            for t in targets:
                if t.chr not in f.references:
                    continue
                for r, al in site_alleles(f, t):
                    if r.mapping_quality < MAPQ:
                        continue
                    kind = "supp" if r.is_supplementary else ("sec" if r.is_secondary else "prim")
                    ub = r.get_tag("UB") if r.has_tag("UB") else None
                    bx = r.get_tag("BX") if r.has_tag("BX") else ub
                    ug = r.get_tag("UG") if r.has_tag("UG") else None
                    rows.append((cell, t.tid, r.query_name, ub, bx, ug, r.reference_start + 1,
                                 r.reference_end, r.is_reverse, kind, al))
        # validation: bam-readcount per-library counts
        cpath = os.path.join(cd, "allcounts.count")
        if os.path.exists(cpath):
            with open(cpath) as fh:
                for line in fh:
                    f_ = line.rstrip("\n").split("\t")
                    chrom, pos, ref = f_[0], int(f_[1]), f_[2]
                    hit = T[(T.chr == chrom) & (T.start == pos) & (T.ref == ref)]
                    if hit.empty:
                        continue
                    for lb, block in re.findall(r"(\S+)\s*\{([^}]*)\}", line):
                        lb = lb.strip("'")
                        for t in hit.itertuples(index=False):
                            c = Counter()
                            for e in block.split():
                                if ":" not in e:
                                    continue
                                b, n = e.split(":")[:2]
                                n = int(n)
                                if b == "=" or b == "N":
                                    continue
                                if b == t.ref:
                                    c["REF"] += n
                                elif t.kind == "SNV" and b == t.alt:
                                    c["ALT"] += n
                                elif t.kind == "indel" and b[0] == t.alt[0] and abs(len(b) - len(t.alt)) <= 1:
                                    c["ALT"] += n
                                else:
                                    c["OTH"] += n
                            rc_groups[(cell, lb, t.tid)] = c

    R = pd.DataFrame(rows, columns=["cell", "tid", "qname", "ub", "bx", "ug", "pos", "end",
                                    "rev", "kind", "allele"])
    R["win100"] = (R.pos // WIN) * WIN
    R["win1k"] = (R.pos // 1000) * 1000
    out = {}

    # ---------------- check 6: secondary / supplementary
    s6 = pd.Series({f"{k}_{'q20' if q else 'lowq'}": v for (k, q), v in rec_stats.items()})
    s6.to_csv(os.path.join(a.out, "c6_record_kinds_all.csv"))
    out["c6_site_records"] = R.groupby(["tid", "kind"]).size().unstack(fill_value=0)

    # ---------------- check 1: bundles
    bs = pd.Series(bundles)
    nu = pd.Series({k: len(v) for k, v in bundle_umis.items()})
    out["c1_bundles"] = pd.DataFrame({
        "bundles": [len(bs)],
        "frac_bundles_1read": [(bs == 1).mean()],
        "frac_bundles_1umi": [(nu == 1).mean()],
        "median_reads_per_bundle": [bs.median()],
        "reads_grouped": [n_grouped],
        "frac_reads_BX_ne_UB": [n_bx_ne_ub / max(n_grouped, 1)],
    })

    # ---------------- groupings
    # current: group key (cell, bx, chr(implicit in tid's chr), win100); size counted over ALL q20
    # records of that key in the cell (as in the awk split), not only those covering the site.
    # The record table only holds site-covering records; window size must count all q20 records
    # in filtered_input, which all overlap some target site -> identical.
    Rd = R.drop_duplicates(["cell", "qname", "kind", "pos", "tid"])
    allrec = Rd.drop_duplicates(["cell", "qname", "kind", "pos"])  # one per record
    tchr = dict(zip(T.tid, T.chr))
    allrec = allrec.assign(chr=allrec.tid.map(tchr))
    cur_size = allrec.groupby(["cell", "bx", "chr", "win100"]).size().rename("n_key")

    def consensus_table(df, keycols, nmin, size=None):
        """Per (cell, tid, group) consensus label; size = Series of group sizes for nmin gate."""
        cnt = df[df.allele.notna()].groupby(["cell", "tid"] + keycols + ["allele"]).size().unstack(fill_value=0)
        for c in ["REF", "ALT", "OTH"]:
            if c not in cnt:
                cnt[c] = 0
        cnt = cnt[["REF", "ALT", "OTH"]]
        if size is None:
            n = df.groupby(["cell", "tid"] + keycols).size().reindex(cnt.index).values
        else:
            n = size.values
        tot = cnt.sum(axis=1).values
        top = cnt.values.argmax(axis=1)
        frac = cnt.values.max(axis=1) / np.maximum(tot, 1)
        lab = np.array(["WT", "MUT", "OTH"])[top]
        ok = (n >= nmin) & (tot > 0) & (frac >= FMIN)
        res = pd.DataFrame({"n": n, "tot": tot, "label": np.where(ok, lab, "drop")}, index=cnt.index)
        return res

    # current grouping (validated against bam-readcount below)
    Rc = Rd.assign(chr=Rd.tid.map(tchr))
    groupings = {}
    cnt_cur = Rc[Rc.allele.notna()].groupby(["cell", "tid", "bx", "chr", "win100", "allele"]).size().unstack(fill_value=0)
    for c in ["REF", "ALT", "OTH"]:
        if c not in cnt_cur:
            cnt_cur[c] = 0
    keys = cnt_cur.index.droplevel("tid")
    nkey = cur_size.reindex(pd.MultiIndex.from_arrays([keys.get_level_values(i) for i in range(4)])).values

    def cur_labels(nmin):
        tot = cnt_cur[["REF", "ALT", "OTH"]].sum(axis=1).values
        v = cnt_cur[["REF", "ALT", "OTH"]].values
        frac = v.max(axis=1) / np.maximum(tot, 1)
        lab = np.array(["WT", "MUT", "OTH"])[v.argmax(axis=1)]
        ok = (nkey >= nmin) & (tot > 0) & (frac >= FMIN)
        return pd.DataFrame({"n": nkey, "tot": tot, "label": np.where(ok, lab, "drop")}, index=cnt_cur.index)

    # ---------------- V: emulation vs bam-readcount (current grouping, nmin=2 gate)
    if rc_groups:
        cl = cnt_cur.copy()
        cl["lb"] = [f"{bx}_{ch}_{w}" for (_, _, bx, ch, w) in cl.index]
        cl["n"] = nkey
        cl = cl[cl.n >= 2]
        agree = tot = 0
        lab_agree = 0
        for (cell, tid, *_), r in cl.iterrows():
            rc = rc_groups.get((cell, r.lb, tid))
            if rc is None:
                continue
            tot += 1
            e = (r.REF, r.ALT, r.OTH)
            b = (rc["REF"], rc["ALT"], rc["OTH"])
            agree += e == b
            le = label({"REF": e[0], "ALT": e[1], "OTH": e[2]}, 1)
            lb_ = label({"REF": b[0], "ALT": b[1], "OTH": b[2]}, 1)
            lab_agree += le == lb_
        out["V_emulation"] = pd.DataFrame({"groups_compared": [tot],
                                           "count_identical": [agree / max(tot, 1)],
                                           "label_identical": [lab_agree / max(tot, 1)],
                                           "rc_groups_total": [len(rc_groups)]})

    # B0: (cell, tid, raw UB), all q20 records; B1: primary only (-F 0x900)
    R1 = Rd[Rd.kind == "prim"]
    # B+C: directional clustering within (cell, tid) on primary records
    def cluster(df, maxd, dist, name):
        reps = {}
        for (cell, tid), sub in df.groupby(["cell", "tid"]):
            u = sub.ub.value_counts().to_dict()
            for k, v in directional(u, maxd, dist).items():
                reps[(cell, tid, k)] = v
        return df.assign(**{name: [reps[(c, t, x)] for c, t, x in zip(df.cell, df.tid, df.ub)]})

    R1 = cluster(R1, 1, hamming, "c_h1")
    R1 = cluster(R1, 1, levenshtein, "c_ld1")
    R1 = cluster(R1, 2, levenshtein, "c_ld2")
    R1 = cluster(R1, 3, levenshtein, "c_ld3")

    schemes = {
        "current(BX x win100, all recs)": lambda n: cur_labels(n),
        "A(current, -F0x900)": lambda n: consensus_table(R1.assign(chr=R1.tid.map(tchr)), ["bx", "chr", "win100"], n),
        "B0(cell x target x UB, all recs)": lambda n: consensus_table(Rd, ["ub"], n),
        "B1(cell x target x UB, -F0x900)": lambda n: consensus_table(R1, ["ub"], n),
        "B1+C(Hamming<=1)": lambda n: consensus_table(R1, ["c_h1"], n),
        "B1+C(LD<=1)": lambda n: consensus_table(R1, ["c_ld1"], n),
        "B1+C(LD<=2)": lambda n: consensus_table(R1, ["c_ld2"], n),
        "B1+C(LD<=3)": lambda n: consensus_table(R1, ["c_ld3"], n),
    }

    # ---------------- check 5: n_min removal and per-target UMI totals / calls
    summ, calls, trans = [], {}, None
    cellcounts = {}
    for name, fn in schemes.items():
        for nmin in NMINS:
            L = fn(nmin)
            labs = L.label.value_counts()
            summ.append({"scheme": name, "nmin": nmin, "groups": len(L),
                         "frac_groups_dropped": (L.label == "drop").mean(),
                         "frac_groups_below_nmin": (L.n < nmin).mean(),
                         "frac_reads_in_groups_below_nmin": L.loc[L.n < nmin, "tot"].sum() / max(L.tot.sum(), 1),
                         "WT": int(labs.get("WT", 0)), "MUT": int(labs.get("MUT", 0)),
                         "OTH": int(labs.get("OTH", 0))})
            tab = L[L.label != "drop"].groupby(["cell", "tid", "label"]).size().unstack(fill_value=0)
            for c in ["WT", "MUT", "OTH"]:
                if c not in tab:
                    tab[c] = 0
            tab = tab[["WT", "MUT", "OTH"]].reset_index()
            # full cell x target grid
            cells = Rd.cell.unique()
            grid = pd.MultiIndex.from_product([cells, T.tid], names=["cell", "tid"]).to_frame(index=False)
            tab = grid.merge(tab, on=["cell", "tid"], how="left").fillna(0)
            ct = call_main(tab)
            cellcounts[(name, nmin)] = ct.set_index(["cell", "tid"])
            calls[(name, nmin)] = ct.groupby(["tid", "genotype"]).size().unstack(fill_value=0)
            pt = ct.groupby("tid")[["WT", "MUT", "OTH"]].sum()
            pt.columns = [f"{c}" for c in pt.columns]
            calls[(name, nmin, "umis")] = pt
    out["c5_nmin_summary"] = pd.DataFrame(summ)
    callrows = []
    for k, v in calls.items():
        if len(k) == 2:
            for tid, r in v.iterrows():
                callrows.append({"scheme": k[0], "nmin": k[1], "tid": tid, **r.to_dict()})
    out["c5_calls_per_target"] = pd.DataFrame(callrows).fillna(0)
    umirows = []
    for k, v in calls.items():
        if len(k) == 3:
            for tid, r in v.iterrows():
                umirows.append({"scheme": k[0], "nmin": k[1], "tid": tid, **r.to_dict()})
    out["c5_umis_per_target"] = pd.DataFrame(umirows)
    # transitions current(nmin=2) -> others
    base = cellcounts[("current(BX x win100, all recs)", 2)]
    tr = []
    for (name, nmin), ct in cellcounts.items():
        j = base[["genotype", "WT", "MUT"]].join(ct[["genotype", "WT", "MUT"]], rsuffix="_new")
        x = pd.crosstab(j.genotype, j.genotype_new)
        for g0 in x.index:
            for g1 in x.columns:
                tr.append({"scheme": name, "nmin": nmin, "from_current_n2": g0, "to": g1, "n": int(x.loc[g0, g1])})
        cov = j[(j.WT + j.MUT) > 0]
        tr.append({"scheme": name, "nmin": nmin, "from_current_n2": "median WT+MUT ratio new/current (covered)",
                   "to": "", "n": float(np.median((cov.WT_new + cov.MUT_new) / (cov.WT + cov.MUT)))})
    out["c5_call_transitions"] = pd.DataFrame(tr)

    # ---------------- check 2: splitting of (cell, tid, UB)
    s2 = []
    for strand_name, sub in [("all", Rd), ("fwd", Rd[~Rd.rev]), ("rev", Rd[Rd.rev])]:
        g = sub.groupby(["cell", "tid", "ub"])
        d = pd.DataFrame({"n": g.size(), "w100": g.win100.nunique(), "w1k": g.win1k.nunique(),
                          "ug": g.ug.nunique(), "bx": g.bx.nunique()})
        for thr in (1, 2, 3):
            dd = d[d.n >= thr]
            s2.append({"strand": strand_name, "min_reads": thr, "ub_molecules": len(dd),
                       "frac_multi_win100": (dd.w100 > 1).mean(), "mean_win100": dd.w100.mean(),
                       "frac_multi_win1k": (dd.w1k > 1).mean(),
                       "frac_multi_UG": (dd.ug > 1).mean(), "mean_UG": dd.ug.mean(),
                       "frac_multi_BX": (dd.bx > 1).mean()})
    out["c2_splitting"] = pd.DataFrame(s2)
    # per target
    g = Rd.groupby(["cell", "tid", "ub"])
    d = pd.DataFrame({"n": g.size(), "w100": g.win100.nunique(), "ug": g.ug.nunique()}).reset_index()
    d = d[d.n >= 2]
    out["c2_splitting_per_target"] = d.groupby("tid").agg(ub_molecules=("n", "size"),
                                                          frac_multi_win100=("w100", lambda x: (x > 1).mean()),
                                                          frac_multi_UG=("ug", lambda x: (x > 1).mean()))
    # read start / end spread within exact-UB molecules (primary, >=3 reads)
    sp = []
    for (cell, tid, ub), sub in R1.groupby(["cell", "tid", "ub"]):
        if len(sub) < 3:
            continue
        five = np.where(sub.rev, sub.end, sub.pos)
        three = np.where(sub.rev, sub.pos, sub.end)
        sp.append((tid, len(sub), np.ptp(sub.pos), np.ptp(sub.end),
                   Counter(sub.pos).most_common(1)[0][1] / len(sub),
                   Counter(sub.end).most_common(1)[0][1] / len(sub),
                   (np.abs(sub.pos - np.median(sub.pos)) <= 5).mean(),
                   (np.abs(sub.end - np.median(sub.end)) <= 5).mean(), sub.rev.mean()))
    sp = pd.DataFrame(sp, columns=["tid", "n", "range_start", "range_end", "frac_modal_start",
                                   "frac_modal_end", "frac_start_within5", "frac_end_within5", "frac_rev"])
    out["c2_end_spread"] = sp.groupby("tid").median().assign(molecules=sp.groupby("tid").size())

    # ---------------- check 3: distance spectrum
    ubs = R1.groupby(["cell", "tid"]).ub.apply(lambda x: sorted(set(x)))
    lens = Counter(len(u) for us in ubs for u in us)
    out["c3_ub_lengths"] = pd.Series(lens).sort_index().to_frame("n")
    spec = []
    for tid in T.tid:
        sub = ubs[ubs.index.get_level_values("tid") == tid]
        within = Counter()
        npairs = 0
        for us in sub:
            npairs += len(us) * (len(us) - 1) // 2
            for i, j, ld in close_pairs(us, 3, levenshtein):
                within["H1"] += hamming(us[i], us[j]) == 1
                within[f"LD{ld}"] += 1
        pool = [(c, u) for (c, _), us in sub.items() for u in us]
        null = Counter()
        nn = 0
        if len(pool) > 1 and sub.size > 1:
            for _ in range(a.null_pairs):
                (c1, u1), (c2, u2) = random.sample(pool, 2)
                if c1 == c2:
                    continue
                nn += 1
                h = hamming(u1, u2)
                ld = levenshtein(u1, u2)
                null["H1"] += h == 1
                null["LD1"] += ld == 1
                null["LD2"] += ld == 2
                null["LD3"] += ld == 3
        row = {"tid": tid, "cells": len(sub), "ubs": len(pool), "within_pairs": npairs, "null_pairs": nn}
        for k in ["H1", "LD1", "LD2", "LD3"]:
            fw = within[k] / max(npairs, 1)
            fn = null[k] / max(nn, 1)
            row[f"{k}_within"] = fw
            row[f"{k}_null"] = fn
            row[f"{k}_excess_pairs"] = within[k] - fn * npairs
        spec.append(row)
    out["c3_distance_spectrum"] = pd.DataFrame(spec)

    # ---------------- check 4: allele concordance of candidate merges
    ublab = consensus_table(R1, ["ub"], 1)
    ublab = ublab[ublab.label.isin(["WT", "MUT"])].reset_index()[["cell", "tid", "ub", "label", "n"]]
    c4 = []
    for tid in T.tid:
        sub = ublab[ublab.tid == tid]
        acc = Counter()
        for cell, cs in sub.groupby("cell"):
            if cs.label.nunique() < 2:
                continue
            ul, ll = list(cs.ub), list(cs.label)
            close = {(i, j): d for i, j, d in close_pairs(ul, 3, levenshtein)}
            for (i, j), d in close.items():
                acc[(f"LD{d}", "pairs")] += 1
                acc[(f"LD{d}", "disc")] += ll[i] != ll[j]
            for _ in range(min(2000, len(ul) * (len(ul) - 1) // 2)):
                i, j = random.sample(range(len(ul)), 2)
                if (min(i, j), max(i, j)) in close:
                    continue
                acc[("LD>3(random)", "pairs")] += 1
                acc[("LD>3(random)", "disc")] += ll[i] != ll[j]
        for cls in ["LD1", "LD2", "LD3", "LD>3(random)"]:
            p = acc[(cls, "pairs")]
            c4.append({"tid": tid, "class": cls, "pairs": p, "discordant": acc[(cls, "disc")],
                       "disc_rate": acc[(cls, "disc")] / p if p else np.nan})
    out["c4_concordance"] = pd.DataFrame(c4)

    for k, v in out.items():
        v.to_csv(os.path.join(a.out, f"{k}.csv"))
    T[["tid", "kind"]].to_csv(os.path.join(a.out, "targets_public.csv"), index=False)
    print("written", sorted(out))


if __name__ == "__main__":
    main()
