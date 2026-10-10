"""Measure the indel labelling change at the KIT deletion target on PALM14E (wayfinder ticket #11).

Decided rule (Indel alleles under the uniform fmin consensus rule, #5): ALT = SUM of indel entries
of the target's type, length within +-k bp, agreeing with the target sequence over their shared
prefix; REF = ref base minus '+' entries; OTH = the rest except N; uniform nmin/fmin, no
shortcuts. main@c442a32 instead takes the MAX over substring matches and calls MUT on any
library with > 3 such reads (len(alt) > 3 shortcut).

bam-readcount -b 20 also filters deletions, on the quality of the anchor base (the read base
before the deletion; checked against bam-readcount -b 20 / -b 0 in brc_debug.sh / brc_delq.py).
The read-level emulation applies the same filter.

Inputs stay on the HPC (cells/<CB>/allcounts.count + filtered_input.bam from stage1, records.pkl
from followup.py, counts_rerun.csv, palm14e_rerun.csv). Writes AGGREGATE tables only; the target
is reported as T3.

  K0 validation: main's rule re-applied to allcounts.count reproduces counts_rerun.csv; read-level
     emulation reproduces bam-readcount per library; main's caller port reproduces its calls
  K1 UMI labels and transitions between rules, main's grouping (bam-readcount libraries) and the
     decided new grouping (cell x target, primary, LD<=2 hybrid T=3, n_min=3)
  K2 cell calls under each rule: main's caller (genotyping.py@c442a32) and the spec caller
     (per-target pooled error, d = 1, floor 1e-4, binomial upper tail, 0.1 / 0.9)
  K3 near-miss spectrum: indel entries on the site line, with and without the anchor-quality filter
  K6 ALT read fraction of UMIs that main's len>3 shortcut called MUT and the new rule calls WT
  K4 blind spot: deletion-like reads (net loss of 3 +- 1 bp in a +-3 bp window) by where the
     aligner put them and what the site line records for them
  window_k1 rule (realignment-lite): ALT = deletion-like read, REF = no indel in the window and REF
     at the site; same quality filters as the site line (site base / deletion anchor >= 20)
"""
import argparse
import glob
import os
import re
from collections import Counter

import numpy as np
import pandas as pd
import pysam
from scipy.stats import binom, nbinom

import hybrid
from measure_umi import FMIN, load_targets

MAPQ = 20
BQ = 20
WIN = 100
EPS = 1e-4
NMIN_MAIN = 2  # baseline params of the PALM14E run (min_read 2, min_fraction 0.6, window 100)
NMIN_NEW = 3   # UMI grouping redesign (#12)
FLANK = 3      # window for deletion-like reads: [site - FLANK, site + len + FLANK)
KS = (0, 1, 2)
LABELS = ["WT", "MUT", "OTH"]


# ---------------------------------------------------------------- rules on entry counters
def main_label(e, ref, alt, nmin=NMIN_MAIN, fmin=FMIN):
    """count_consensus.consensus_on_umi of main@c442a32 on one library's entries."""
    c = {"ref": 0, "alt": 0, "mis": 0}
    for b, n in e.items():
        if b == ref:
            c["ref"] = n
        elif len(alt) <= 2 and b == alt:
            c["alt"] = n
        elif len(alt) > 2 and len(b) > 2 and (b in alt or alt in b):
            c["alt"] = max(n, c["alt"])
        elif b != "N":
            c["mis"] += n
    tot = sum(c.values())
    if len(alt) == 2 and c["alt"] > 2 and c["alt"] / (c["ref"] + c["alt"]) > 0.25:
        return "MUT"
    if len(alt) > 3 and c["alt"] > 3:
        return "MUT"
    m = max(c["ref"], c["alt"])
    if m >= nmin and m / max(tot, 1) >= fmin:
        return {"ref": "WT", "alt": "MUT", "mis": "OTH"}[max(c, key=c.get)]
    return "drop"


def shortcut_mut(e, alt):
    return main_label(e, "\0", alt, nmin=10**9) == "MUT"


def matches(b, alt, k):
    """Tolerant indel match: same type, |len diff| <= k, agree over the shared prefix."""
    if b[0] != alt[0] or abs(len(b) - len(alt)) > k:
        return False
    s, t = b[1:], alt[1:]
    m = min(len(s), len(t))
    return s[:m] == t[:m]


def new_counts(e, ref, alt, k):
    plus = sum(n for b, n in e.items() if b.startswith("+"))
    c = Counter()
    for b, n in e.items():
        if b in ("=", "N") or n == 0:
            continue
        if b == ref:
            c["REF"] += n
        elif b[0] in "+-" and matches(b, alt, k):
            c["ALT"] += n
        else:
            c["OTH"] += n
    c["REF"] -= min(plus, c["REF"])  # reads with an insertion are also in the anchor's base count
    return c


def uniform_label(c, n_group, nmin, top_min):
    """Uniform rule: group gate n_group >= nmin; top allele >= top_min reads and >= fmin."""
    tot = c["REF"] + c["ALT"] + c["OTH"]
    if n_group < nmin or tot == 0:
        return "drop"
    top = max(("REF", "ALT", "OTH"), key=lambda x: c[x])
    if c[top] < top_min or c[top] / tot < FMIN:
        return "drop"
    return {"REF": "WT", "ALT": "MUT", "OTH": "OTH"}[top]


# ---------------------------------------------------------------- callers
def call_main_exact(df):
    """genotyping.py@c442a32: per-cell err = OTH / all / 2 over targets; err == 0 -> min(err) (0 when
    any cell has 0); floor 1e-4; nbinom.cdf(WT, MUT, err); MUT < 0.1, WT >= 0.9 or (NaN, WT >= 2)."""
    g = df.groupby("cell")
    err = (g["OTH"].sum() / (g[LABELS].sum().sum(axis=1) + 1e-8) / 2).rename("err")
    df = df.merge(err, left_on="cell", right_index=True)
    df.loc[df.err == 0, "err"] = df.err.min()
    e = np.maximum(df.err.values, EPS)
    wt, mut = df.WT.values, df.MUT.values
    with np.errstate(all="ignore"):
        p = nbinom.cdf(k=wt, n=mut, p=e)
    gt = np.full(p.size, "NA", dtype=object)
    gt[p < 0.1] = "MUT"
    gt[np.isnan(p) & (wt >= 2)] = "WT"
    gt[p >= 0.9] = "WT"
    gt[(wt == 0) & (mut == 0)] = "NA"
    return df.assign(genotype=gt)


def call_spec(tab):
    """Revised document + #6 for one indel target: P_s = max(w / N, eps) pooled over cells
    (d = 1); p = P(X >= m | Binom(m + k, P_s)); MUT m >= 1 & p < 0.1; WT (m = 0 & k >= 2) or p >= 0.9."""
    w, N = tab.OTH.sum(), tab[LABELS].values.sum()
    ps = max(w / N if N else 0, EPS)
    m, k = tab.MUT.values, tab.WT.values
    p = binom.sf(m - 1, m + k, ps)
    gt = np.full(len(tab), "NA", dtype=object)
    gt[(m >= 1) & (p < 0.1)] = "MUT"
    gt[((m == 0) & (k >= 2)) | ((m >= 1) & (p >= 0.9))] = "WT"
    return tab.assign(genotype=gt), ps


# ---------------------------------------------------------------- read level
def read_info(r, t, fa_seq, w0):
    """Site-line entries bam-readcount would record for this read (raw and after its quality
    filter), plus the window class. fa_seq: reference over [w0, ...) (0-based start w0)."""
    s0 = t.start - 1
    L = len(t.alt) - 1
    lo, hi = s0 - FLANK, s0 + L + FLANK
    q = r.query_qualities
    dels, ins = [], []  # (ref start, len, anchor quality) / (anchor ref pos, seq)
    site = None
    rp, qp = r.reference_start, 0
    for op, ln in r.cigartuples:
        if op in (0, 7, 8):  # M = X
            if rp <= s0 < rp + ln:
                site = ("base", qp + s0 - rp)
            rp += ln
            qp += ln
        elif op == 1:  # I, anchored at rp - 1
            ins.append((rp - 1, r.query_sequence[qp:qp + ln]))
            qp += ln
        elif op in (2, 3):  # D / N
            aq = q[qp - 1] if qp > 0 else 0
            if op == 2:
                dels.append((rp, ln, aq))
                if rp == s0:
                    site = ("del", ln, aq)
            if site is None and rp <= s0 < rp + ln:
                site = ("in_del",)
            rp += ln
        elif op == 4:
            qp += ln
    raw, brc = (), ()
    if site is not None and site[0] == "base":
        e = (r.query_sequence[site[1]],) + tuple("+" + s for a, s in ins if a == s0)
        raw = e
        if q[site[1]] >= BQ:
            brc = e
    elif site is not None and site[0] == "del":
        raw = ("-" + fa_seq[s0 - w0:s0 - w0 + site[1]],)
        if site[2] >= BQ:
            brc = raw
    covers = r.reference_start <= lo and r.reference_end >= hi
    wd = [d for d in dels if d[0] < hi and d[0] + d[1] > lo]
    del_bp = sum(min(a + b, hi) - max(a, lo) for a, b, _ in wd)
    ins_bp = sum(len(s) for a, s in ins if lo <= a < hi)
    if not wd:
        wclass = "no_del"
    elif len(wd) > 1:
        wclass = "split"
    elif wd[0][:2] == (s0, L):
        wclass = "exact_site"
    elif wd[0][0] == s0:
        wclass = "site_other_len"
    else:
        wclass = f"shifted_{wd[0][0] - s0:+d}"
    wq = wd[0][2] if wd else -1
    return raw, brc, covers, wclass, del_bp - ins_bp, del_bp, ins_bp, wq


def site_verdict(entries, t, k=1):
    """What the site line makes of this read under the new rule."""
    if not entries:
        return "absent"
    b = entries[0]
    if b[0] == "-":
        return "ALT" if matches(b, t.alt, k) else "OTH_indel"
    if len(entries) > 1:
        return "OTH_ins"
    return "REF" if b == t.ref else ("absent" if b == "N" else "OTH_base")


def window_allele(row, t, k=1):
    """Realignment-lite with the site line's quality filters (see module doc); reads not spanning
    the window fall back to the site line."""
    if not row.covers:
        v = site_verdict(row.brc, t, k)
        return {"ALT": "ALT", "REF": "REF", "absent": None}.get(v, "OTH")
    L = len(t.alt) - 1
    if row.del_bp > 0 and abs(row.net - L) <= k:
        return "ALT" if row.wq >= BQ else None
    if row.del_bp == 0 and row.ins_bp == 0:
        if not row.brc:
            return None
        return "REF" if row.brc[0] == t.ref else ("OTH" if row.brc[0] != "N" else None)
    return "OTH" if (row.del_bp == 0 or row.wq >= BQ) else None


def parse_line(line):
    """bam-readcount -p line -> {library: Counter(entry -> reads)}."""
    out = {}
    for lb, block in re.findall(r"(\S+)\s*\{([^}]*)\}", line):
        c = Counter()
        for e in block.split():
            if ":" in e:
                b, n = e.split(":")[:2]
                if int(n):
                    c[b] += int(n)
        out[lb.strip("'")] = c
    return out


def label_table(index_cells, labels, cells, tid):
    """Per-cell WT/MUT/OTH UMI counts at the target from per-group labels."""
    s = pd.DataFrame({"cell": np.asarray(index_cells), "label": np.asarray(labels)})
    s = s[s.label != "drop"]
    tab = s.groupby(["cell", "label"]).size().unstack(fill_value=0).reindex(columns=LABELS, fill_value=0)
    return tab.reindex(cells, fill_value=0).rename_axis("cell").reset_index().assign(tid=tid)


def calls_for(tab, others, tid):
    cm = call_main_exact(pd.concat([others, tab], ignore_index=True))
    cm = cm[cm.tid == tid].set_index("cell")
    cs, ps = call_spec(tab.set_index("cell"))
    return {"main_caller": cm, "spec_caller": cs}, ps


def summarise_calls(grouping, rule, calls, ps, k2):
    for caller, c in calls.items():
        cov = c[LABELS].sum(axis=1) > 0
        k2.append({"grouping": grouping, "rule": rule, "caller": caller,
                   "cells_with_umis": int(cov.sum()),
                   "MUT": int((c.genotype == "MUT").sum()), "WT": int((c.genotype == "WT").sum()),
                   "NA_with_umis": int(((c.genotype == "NA") & cov).sum()),
                   "MUT_from_single_MUT_umi": int(((c.genotype == "MUT") & (c.MUT == 1)).sum()),
                   "umis_WT": int(c.WT.sum()), "umis_MUT": int(c.MUT.sum()), "umis_OTH": int(c.OTH.sum()),
                   "P_s": ps if caller == "spec_caller" else np.nan})


def transitions(grouping, pairs, labels, calls, k1, k2t):
    for fr, to in pairs:
        for (f, x), n in pd.crosstab(labels[fr], labels[to]).stack().items():
            if f != x and n:
                k1.append({"grouping": grouping, "from_rule": fr, "to_rule": to, "from": f, "to": x, "umis": int(n)})
        for caller in calls[fr]:
            j = calls[fr][caller][["genotype"]].join(calls[to][caller][["genotype"]], rsuffix="_to")
            for (f, x), n in j.groupby(["genotype", "genotype_to"]).size().items():
                if f != x:
                    k2t.append({"grouping": grouping, "caller": caller, "from_rule": fr, "to_rule": to,
                                "from": f, "to": x, "cells": int(n)})


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--counts", required=True, help="counts_rerun.csv (main's per-cell table)")
    ap.add_argument("--calls", required=True, help="palm14e_rerun.csv (main's calls on the re-run)")
    ap.add_argument("--cache", default="records.pkl")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    T = load_targets(a.bed)
    t = next(T[T.kind == "indel"].itertuples(index=False))
    tid = t.tid
    w0 = t.start - 1 - 50
    fa_seq = pysam.FastaFile(a.ref).fetch(t.chr, w0, t.start - 1 + 50).upper()

    libs = {}   # (cell, lib) -> Counter from bam-readcount
    reads = []  # one per q20 record overlapping the window
    for cd in sorted(glob.glob(os.path.join(a.cells, "*"))):
        cell = os.path.basename(cd)
        cp = os.path.join(cd, "allcounts.count")
        if os.path.exists(cp):
            with open(cp) as fh:
                for line in fh:
                    f = line.split("\t", 3)
                    if f[0] == t.chr and int(f[1]) == t.start and f[2] == t.ref:
                        for lb, c in parse_line(line).items():
                            libs[(cell, lb)] = c
        bam = os.path.join(cd, "filtered_input.bam")
        if not os.path.exists(bam + ".bai"):
            bam = os.path.join(cd, "sorted_fi.bam")
        if not os.path.exists(bam):
            continue
        with pysam.AlignmentFile(bam, "rb") as fb:
            if t.chr not in fb.references:
                continue
            for r in fb.fetch(t.chr, t.start - 1 - FLANK, t.start + FLANK + 3):
                if r.is_unmapped or r.mapping_quality < MAPQ or r.cigartuples is None:
                    continue
                kind = "supp" if r.is_supplementary else ("sec" if r.is_secondary else "prim")
                bx = r.get_tag("BX") if r.has_tag("BX") else None
                reads.append((cell, r.query_name, r.reference_start + 1, kind,
                              f"{bx}_{t.chr}_{(r.reference_start + 1) // WIN * WIN}",
                              *read_info(r, t, fa_seq, w0)))
    Rd = pd.DataFrame(reads, columns=["cell", "qname", "pos", "kind", "lib", "raw", "brc", "covers",
                                      "wclass", "net", "del_bp", "ins_bp", "wq"])
    Rd["site_raw"] = [site_verdict(e, t, 1) for e in Rd.raw]
    Rd["site_brc"] = [site_verdict(e, t, 1) for e in Rd.brc]
    Rd["win_allele"] = [window_allele(r, t, 1) for r in Rd.itertuples(index=False)]
    print("reads", len(Rd), "libraries", len(libs), flush=True)

    # ---------------- main's grouping: bam-readcount libraries
    L = pd.DataFrame([{"cell": c, "lib": lb, "entries": e} for (c, lb), e in libs.items()])
    cr = pd.read_csv(a.counts)
    tmap = {(r.chr, r.start, r.alt): r.tid for r in T.itertuples()}
    cr["tid"] = [tmap[(c, s, x)] for c, s, x in zip(cr.chr, cr.start, cr.alt)]
    cr = cr.rename(columns={"barcode": "cell", "MIS": "OTH"})
    cells = cr.cell.unique()
    others = cr[cr.tid != tid][["cell", "tid"] + LABELS]

    L["main"] = [main_label(e, t.ref, t.alt) for e in L.entries]
    for k in KS:
        L[f"new_k{k}"] = [uniform_label(new_counts(e, t.ref, t.alt, k), NMIN_MAIN, NMIN_MAIN, NMIN_MAIN)
                          for e in L.entries]
    key = list(zip(L.cell, L.lib))
    emu = Rd.groupby(["cell", "lib"]).brc.agg(lambda s: Counter(x for e in s for x in e))
    emu_raw = Rd.groupby(["cell", "lib"]).raw.agg(lambda s: Counter(x for e in s for x in e))
    wl = Rd.groupby(["cell", "lib"]).win_allele.agg(lambda s: Counter(x for x in s if isinstance(x, str)))
    L["site_k1_noqual_del"] = [uniform_label(new_counts(emu_raw.get(x, Counter()), t.ref, t.alt, 1),
                                             NMIN_MAIN, NMIN_MAIN, NMIN_MAIN) for x in key]
    L["window_k1"] = [uniform_label(wl.get(x, Counter()), NMIN_MAIN, NMIN_MAIN, NMIN_MAIN) for x in key]
    L["shortcut"] = [shortcut_mut(e, t.alt) for e in L.entries]

    # K0 validation
    main_tab = label_table(L.cell, L["main"], cells, tid).set_index("cell")[LABELS]
    crt = cr[cr.tid == tid].set_index("cell")[LABELS].reindex(main_tab.index)
    ident = [emu.get(x, Counter()) == Counter({b: n for b, n in e.items() if b != "="})
             for x, e in zip(key, L.entries)]
    base = pd.read_csv(a.calls)
    base = base[(base.start == t.start) & (base.alt == t.alt)].set_index("cell").genotype.fillna("NA")
    port = call_main_exact(cr[["cell", "tid"] + LABELS].copy())
    port = port[port.tid == tid].set_index("cell").genotype
    pd.DataFrame({
        "cells_with_site_line": [L.cell.nunique()], "libraries": [len(L)],
        "cells_main_rule_reproduces_counts_rerun": [(main_tab == crt).all(axis=1).mean()],
        "libraries_read_emulation_identical": [np.mean(ident)],
        "main_caller_port_reproduces_calls": [(port.reindex(base.index) == base).mean()],
    }).to_csv(os.path.join(a.out, "K0_validation.csv"), index=False)

    rules = ["main"] + [f"new_k{k}" for k in KS] + ["site_k1_noqual_del", "window_k1"]
    pairs = [("main", "new_k0"), ("main", "new_k1"), ("main", "new_k2"),
             ("new_k1", "site_k1_noqual_del"), ("new_k1", "window_k1")]
    k1s, k1, k2, k2t = [], [], [], []
    calls_m = {}
    for r in rules:
        k1s.append({"grouping": "main", "rule": r, **L[r].value_counts().to_dict()})
        calls_m[r], ps = calls_for(label_table(L.cell, L[r], cells, tid), others, tid)
        summarise_calls("main", r, calls_m[r], ps, k2)
    k1s.append({"grouping": "main", "rule": "main: MUT via len>3 shortcut",
                "MUT": int((L.shortcut & (L["main"] == "MUT")).sum())})
    transitions("main", pairs, L, calls_m, k1, k2t)

    # ---------------- K3 near-miss spectrum
    prim = Rd[Rd.kind == "prim"].drop_duplicates(["cell", "qname", "pos"])
    k3 = []
    lib_cnt = Counter()
    for e in L.entries:
        lib_cnt.update({b: n for b, n in e.items() if b != "="})
    for src, cnt in [("bam-readcount libraries", lib_cnt),
                     ("primary reads, site line (anchor q>=20)", Counter(x for e in prim.brc for x in e)),
                     ("primary reads, no quality filter", Counter(x for e in prim.raw for x in e))]:
        tot = sum(cnt.values())
        for b, n in cnt.most_common():
            if b[0] in "+-" or n >= 10:
                k3.append({"source": src, "entry": b, "reads": n, "frac_of_site_entries": n / tot,
                           "main_alt": b[0] in "+-" and len(b) > 2 and (b in t.alt or t.alt in b),
                           **{f"alt_k{k}": b[0] in "+-" and matches(b, t.alt, k) for k in KS}})
    pd.DataFrame(k3).groupby("source").head(25).to_csv(os.path.join(a.out, "K3_near_miss_spectrum.csv"), index=False)

    # ---------------- K4 blind spot (primary reads spanning the window)
    pc = prim[prim.covers].copy()
    pc["del_like_k1"] = (pc.del_bp > 0) & ((pc.net - (len(t.alt) - 1)).abs() <= 1)
    pc["wgroup"] = np.where(pc.wclass.str.startswith("shifted"), "shifted", pc.wclass)
    pc.groupby(["del_like_k1", "wgroup", "site_raw"]).size().rename("reads").reset_index() \
      .to_csv(os.path.join(a.out, "K4_blind_spot.csv"), index=False)
    pc[pc.wclass.str.startswith("shifted") & pc.del_like_k1].wclass.value_counts().rename("reads") \
      .rename_axis("shift").reset_index().to_csv(os.path.join(a.out, "K4_shift_offsets.csv"), index=False)
    dl = pc[pc.del_like_k1]
    dq = dl[dl.wq >= BQ]
    pd.DataFrame({
        "primary_reads_spanning_window": [len(pc)], "deletion_like_k1": [len(dl)],
        "frac_on_site_line_ALT_k1": [(dl.site_raw == "ALT").mean()],
        "frac_off_line_absent": [(dl.site_raw == "absent").mean()],
        "frac_off_line_REF": [(dl.site_raw == "REF").mean()],
        "frac_off_line_OTH": [dl.site_raw.str.startswith("OTH").mean()],
        "deletion_like_anchor_q20": [len(dq)],
        "anchor_q20_frac_on_site_line_ALT_k1": [(dq.site_raw == "ALT").mean()],
        "frac_deletion_like_failing_anchor_q20": [1 - len(dq) / max(len(dl), 1)],
        "frac_REF_site_reads_failing_base_q20": [((pc.site_raw == "REF") & (pc.site_brc == "absent")).sum()
                                                 / max((pc.site_raw == "REF").sum(), 1)],
        "non_deletion_like_ALT_k1": [int(((~pc.del_like_k1) & (pc.site_raw == "ALT")).sum())],
    }).to_csv(os.path.join(a.out, "K4_summary.csv"), index=False)

    # ---------------- new grouping: hybrid T=3 clusters of primary records (records.pkl)
    R = pd.read_pickle(a.cache)
    R1 = R[R.kind == "prim"].drop_duplicates(["cell", "qname", "pos", "tid"]).reset_index(drop=True)
    hybrid.SCHEMES = {"T3": (3, True)}
    R1 = hybrid.assign(R1)
    ot = hybrid.consensus(R1[R1.tid != tid], "cl_T3")
    ot = ot[ot.label != "drop"].groupby(["cell", "tid", "label"]).size() \
        .unstack(fill_value=0).reindex(columns=LABELS, fill_value=0).reset_index()
    cells_n = R1.cell.unique()
    grid = pd.MultiIndex.from_product([cells_n, T.tid], names=["cell", "tid"]).to_frame(index=False)
    others_n = grid[grid.tid != tid].merge(ot, on=["cell", "tid"], how="left").fillna(0)
    Rt = R1[R1.tid == tid].merge(prim[["cell", "qname", "pos", "raw", "brc", "win_allele"]],
                                 on=["cell", "qname", "pos"], how="left")
    for c in ("raw", "brc"):
        Rt[c] = [e if isinstance(e, tuple) else () for e in Rt[c]]
    g = Rt.groupby(["cell", "cl_T3"])
    G = pd.DataFrame({"n": g.size(),
                      "brc": g.brc.agg(lambda s: Counter(x for e in s for x in e)),
                      "raw": g.raw.agg(lambda s: Counter(x for e in s for x in e)),
                      "win": g.win_allele.agg(lambda s: Counter(x for x in s if isinstance(x, str)))})
    G["main"] = [main_label(e, t.ref, t.alt, nmin=1) if n >= NMIN_NEW else "drop" for e, n in zip(G.brc, G.n)]
    for k in KS:
        G[f"new_k{k}"] = [uniform_label(new_counts(e, t.ref, t.alt, k), n, NMIN_NEW, 1) for e, n in zip(G.brc, G.n)]
    G["site_k1_noqual_del"] = [uniform_label(new_counts(e, t.ref, t.alt, 1), n, NMIN_NEW, 1) for e, n in zip(G.raw, G.n)]
    G["window_k1"] = [uniform_label(c, n, NMIN_NEW, 1) for c, n in zip(G.win, G.n)]
    G["shortcut"] = [shortcut_mut(e, t.alt) for e in G.brc]
    G = G.reset_index()
    calls_n = {}
    for r in rules:
        k1s.append({"grouping": "new", "rule": r, **G[r].value_counts().to_dict()})
        calls_n[r], ps = calls_for(label_table(G.cell, G[r], cells_n, tid), others_n, tid)
        summarise_calls("new", r, calls_n[r], ps, k2)
    k1s.append({"grouping": "new", "rule": "main: MUT via len>3 shortcut",
                "MUT": int((G.shortcut & (G["main"] == "MUT")).sum())})
    transitions("new", pairs, G, calls_n, k1, k2t)

    # ---------------- K6 ALT read fraction of UMIs flipped MUT -> WT by dropping the shortcut
    k6 = []
    for grouping, D_, col in [("main", L, "entries"), ("new", G, "brc")]:
        for name, m in [("main MUT -> new_k1 WT", (D_["main"] == "MUT") & (D_["new_k1"] == "WT")),
                        ("MUT under both", (D_["main"] == "MUT") & (D_["new_k1"] == "MUT"))]:
            c = [new_counts(e, t.ref, t.alt, 1) for e in D_.loc[m, col]]
            fr = pd.Series([x["ALT"] / max(sum(x.values()), 1) for x in c])
            n = pd.Series([sum(x.values()) for x in c])
            k6.append({"grouping": grouping, "umis": name, "n": int(m.sum()),
                       "reads_median": n.median(), "alt_frac_median": fr.median(),
                       **{f"alt_frac_{lo}-{hi}": int(((fr >= lo) & (fr < hi)).sum())
                          for lo, hi in [(0, .1), (.1, .2), (.2, .3), (.3, .4), (.4, 1.01)]}})
    pd.DataFrame(k6).to_csv(os.path.join(a.out, "K6_shortcut_flips.csv"), index=False)

    pd.DataFrame(k1s).fillna(0).to_csv(os.path.join(a.out, "K1_umi_labels.csv"), index=False)
    pd.DataFrame(k1).to_csv(os.path.join(a.out, "K1_umi_transitions.csv"), index=False)
    pd.DataFrame(k2).to_csv(os.path.join(a.out, "K2_calls.csv"), index=False)
    pd.DataFrame(k2t).to_csv(os.path.join(a.out, "K2_call_transitions.csv"), index=False)
    print("done")


if __name__ == "__main__":
    main()
