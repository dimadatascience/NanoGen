"""Run the main -> v2 attribution ladder on PALM14E (wayfinder ticket #21).

Design: Regression check: main vs genotype-v2 calls on PALM14E (#19), part (b). Cumulative rungs:
  S0  main emulated: "current" grouping (BX x chr x 100 bp window, all q20 records), n_min 2,
      main's consensus rule (count_consensus.py@c442a32), main's caller (call_main_exact)
  S1  + primary records only (-F 0x904)
  S2  + grouping: cell x target, Levenshtein <= 2 hybrid T = 3 (hybrid.py)
  S3  + n_min 2 -> 3
  S4  + spec consensus rule: uniform n_min / fmin with OTH tested on o; indel ALT tolerant k = 1,
      no shortcuts (kit_indel.new_counts / uniform_label). The n_min gate counts reads with a usable
      allele at the site (r + a + o, revised document section 1)
  S4n off-ladder: S4 with the gate on all cluster reads instead (the convention of hybrid.py /
      kit_indel.py), to size that ambiguity
  S5  + indel target at -b 0 (no base / anchor quality filter at the indel site line)
  S6  + spec caller: per-target pooled P_s = max(w / (d N), 1e-4), d = 2 SNV / 1 indel,
      binomial upper tail, 0.1 / 0.9  (= v2 expected)
  S0+6  off-ladder: the spec caller on S0's counts.

Inputs stay on the HPC: records.pkl (followup.py), cells/<CB>/filtered_input.bam for the indel
site line (kit_indel.read_info), counts_rerun.csv / palm14e_rerun.csv / palm14e.csv for the S0 tie.
Writes AGGREGATE tables (targets as T1..Tn) to <out>/aggregate and per-cell / per-flip tables to
<out>/private (HPC only, never committed):
  A0 S0 tie to main's real output        A1 per-rung summary and verdict per target
  A2 transition matrices                  A3 (k, m, o) of hard flips
  A4 per-cell S0 -> S6                    A5 error table (main per-cell err vs v2 pooled P_s)
"""
import argparse
import glob
import os
from collections import Counter

import numpy as np
import pandas as pd
import pysam
from scipy.stats import binom

import hybrid
from kit_indel import FLANK, LABELS, call_main_exact, main_label, new_counts, read_info, uniform_label
from measure_umi import FMIN, MAPQ, WIN, load_targets

EPS = 1e-4
NOISE = 0.001     # acceptance 1: hard flips against prediction <= 0.1 % of covered rows
LOCUS_BP = 1000   # targets closer than this on one chromosome share reads: one independent locus
HARD = ("MUT", "WT")

# rung -> (records, grouping, n_min, consensus rule, indel site line, caller)
RUNGS = {
    "S0": ("all", "main", 2, "main", "b20", "main"),
    "S1": ("prim", "main", 2, "main", "b20", "main"),
    "S2": ("prim", "hybrid", 2, "main", "b20", "main"),
    "S3": ("prim", "hybrid", 3, "main", "b20", "main"),
    "S4": ("prim", "hybrid", 3, "spec", "b20", "main"),
    "S5": ("prim", "hybrid", 3, "spec", "b0", "main"),
    "S6": ("prim", "hybrid", 3, "spec", "b0", "spec"),
    "S0+6": ("all", "main", 2, "main", "b20", "spec"),
    "S4n": ("prim", "hybrid", 3, "spec_n", "b20", "main"),
}
PAIRS = [("S0", "S1"), ("S1", "S2"), ("S2", "S3"), ("S3", "S4"), ("S4", "S5"), ("S5", "S6"),
         ("S0", "S0+6"), ("S0+6", "S6"), ("S0", "S6"), ("S3", "S4n"), ("S4", "S4n")]


def allowed(rung, kind, d):
    """Predicted hard-flip directions (#19 resolution). d is 'MUT>WT' or 'WT>MUT'."""
    if rung in ("S6", "S0+6"):
        return True                      # caller: either direction
    if rung == "S2":
        return d == "MUT>WT"             # fragmented WT molecules
    if rung in ("S4", "S4n"):
        return kind == "indel" and d == "MUT>WT"
    if rung == "S5":
        return kind == "indel"           # few, either direction, indel target only
    return False                         # S1 few, S3 none


# ---------------------------------------------------------------- inputs
def indel_entries(cells, cellset, t, ref):
    """Site-line entries per record at the indel target: b20 = bam-readcount -b 20 (site base and
    deletion anchor >= 20), b0 = no quality filter."""
    w0 = t.start - 1 - 50
    fa_seq = pysam.FastaFile(ref).fetch(t.chr, w0, t.start - 1 + 50).upper()
    rows = []
    for cd in sorted(glob.glob(os.path.join(cells, "*"))):
        cell = os.path.basename(cd)
        if cell not in cellset:
            continue
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
                raw, brc = read_info(r, t, fa_seq, w0)[:2]
                rows.append((cell, r.query_name, r.reference_start + 1, kind, brc, raw))
    E = pd.DataFrame(rows, columns=["cell", "qname", "pos", "kind", "b20", "b0"])
    return E.drop_duplicates(["cell", "qname", "pos", "kind"]).assign(tid=t.tid)


# ---------------------------------------------------------------- counts per rung
def group_labels(D, keycols, nmin, rule, t3, site):
    """Per (cell, tid, group) consensus label. D has column n (group read count for the gate)."""
    keys = ["cell", "tid"] + keycols
    out = []
    S = D[D.tid != t3.tid]
    cnt = S[S.allele.notna()].groupby(keys + ["allele"]).size().unstack(fill_value=0) \
        .reindex(columns=["REF", "ALT", "OTH"], fill_value=0)
    n = S.groupby(keys).n.first().reindex(cnt.index).values
    v = cnt.values
    tot = v.sum(axis=1)
    lab = np.array(LABELS)[v.argmax(axis=1)]
    if rule == "main":   # gate on max(ref, alt) reads; OTH never wins at an SNV
        top = np.maximum(v[:, 0], v[:, 1])
        ok = (n >= nmin) & (top >= nmin) & (top / np.maximum(tot, 1) >= FMIN)
    else:                # uniform: gate on usable reads (spec) or cluster reads (spec_n); top (OTH included) >= fmin
        gate = n if rule == "spec_n" else tot
        ok = (gate >= nmin) & (tot > 0) & (v.max(axis=1) / np.maximum(tot, 1) >= FMIN)
    out.append(pd.DataFrame({"label": np.where(ok, lab, "drop")}, index=cnt.index).reset_index()[["cell", "tid", "label"]])
    I = D[D.tid == t3.tid]
    if len(I):
        g = I.groupby(keys)
        G = pd.DataFrame({"n": g.n.first(), "e": g[site].agg(lambda s: Counter(x for e in s for x in e))})
        if rule == "main":
            G["label"] = [main_label(e, t3.ref, t3.alt, nmin) if n >= nmin else "drop" for e, n in zip(G.e, G.n)]
        else:
            C = [new_counts(e, t3.ref, t3.alt, 1) for e in G.e]
            gates = G.n if rule == "spec_n" else [c["REF"] + c["ALT"] + c["OTH"] for c in C]
            G["label"] = [uniform_label(c, n, nmin, 1) for c, n in zip(C, gates)]
        out.append(G.reset_index()[["cell", "tid", "label"]])
    return pd.concat(out, ignore_index=True)


def cell_table(L, grid):
    tab = L[L.label != "drop"].groupby(["cell", "tid", "label"]).size().unstack(fill_value=0) \
        .reindex(columns=LABELS, fill_value=0).reset_index()
    return grid.merge(tab, on=["cell", "tid"], how="left").fillna(0).astype({c: int for c in LABELS})


def call_spec_all(tab, T):
    """Spec caller per target: P_s pooled over cells, d = 2 (SNV) / 1 (indel), floor 1e-4."""
    out, ps = [], {}
    for t in T.itertuples(index=False):
        s = tab[tab.tid == t.tid]
        w, N = s.OTH.sum(), s[LABELS].values.sum()
        d = 1 if t.kind == "indel" else 2
        q = w / (d * N) if N else 0.0
        ps[t.tid] = (q, max(q, EPS))
        m, k = s.MUT.values, s.WT.values
        p = binom.sf(m - 1, m + k, ps[t.tid][1])
        gt = np.full(len(s), "NA", dtype=object)
        gt[(m >= 1) & (p < 0.1)] = "MUT"
        gt[((m == 0) & (k >= 2)) | ((m >= 1) & (p >= 0.9))] = "WT"
        out.append(s.assign(p=p, genotype=gt))
    return pd.concat(out, ignore_index=True), ps


# ---------------------------------------------------------------- comparisons
def q(x):
    x = pd.Series(x)
    if x.empty:
        return ""
    return "{:g} [{:g}-{:g}] ({:g}-{:g})".format(x.median(), x.quantile(.25), x.quantile(.75), x.min(), x.max())


def compare(a, b, A, B, T, a1, a2, a3, flips):
    j = A.set_index(["cell", "tid"])[LABELS + ["genotype"]].join(
        B.set_index(["cell", "tid"])[LABELS + ["genotype"]], rsuffix="_b")
    for t in T.itertuples(index=False):
        s = j[j.index.get_level_values("tid") == t.tid]
        cov = (s[LABELS].sum(axis=1) > 0) | (s[[c + "_b" for c in LABELS]].sum(axis=1) > 0)
        s = s[cov]
        for (g0, g1), n in s.groupby(["genotype", "genotype_b"]).size().items():
            a2.append({"from_rung": a, "to_rung": b, "tid": t.tid, "from": g0, "to": g1, "rows": int(n)})
        hard = s.genotype.isin(HARD) & s.genotype_b.isin(HARD) & (s.genotype != s.genotype_b)
        soft = (s.genotype != s.genotype_b) & ~hard
        row = {"from_rung": a, "to_rung": b, "tid": t.tid, "kind": t.kind, "covered_rows": int(len(s)),
               "soft_flips": int(soft.sum()), "soft_to_NA": int((soft & (s.genotype_b == "NA")).sum())}
        against = 0
        for d, (g0, g1) in {"MUT>WT": ("MUT", "WT"), "WT>MUT": ("WT", "MUT")}.items():
            h = s[(s.genotype == g0) & (s.genotype_b == g1)]
            row[d] = len(h)
            if not allowed(b, t.kind, d):
                against += len(h)
            unchanged = (h.WT == h.WT_b) & (h.MUT == h.MUT_b) & (h.OTH == h.OTH_b)
            row[f"{d}_counts_unchanged"] = int(unchanged.sum())
            if len(h):
                a3.append({"from_rung": a, "to_rung": b, "tid": t.tid, "direction": d, "flips": len(h),
                           "counts_unchanged": int(unchanged.sum()),
                           **{f"{c}_before": q(h[c]) for c in LABELS},
                           **{f"{c}_after": q(h[c + "_b"]) for c in LABELS}})
                flips.append(h.reset_index().assign(from_rung=a, to_rung=b, direction=d))
        row["against_prediction"] = against
        row["against_frac"] = against / max(len(s), 1)
        if (a, b) in PAIRS[:6] or b in ("S0+6", "S4n"):
            row["verdict_direction"] = "ok" if row["against_frac"] <= NOISE else "BLOCK"
        a1.append(row)


def loci(T):
    """Independent loci: targets within LOCUS_BP on one chromosome share reads."""
    lab, cur, prev = {}, 0, None
    for t in T.sort_values(["chr", "start"]).itertuples(index=False):
        if prev is None or t.chr != prev.chr or t.start - prev.start > LOCUS_BP:
            cur += 1
        lab[t.tid] = f"L{cur}"
        prev = t
    return lab


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--cache", required=True, help="records.pkl")
    ap.add_argument("--counts", required=True, help="counts_rerun.csv (main's counts on this subset)")
    ap.add_argument("--calls", nargs="+", required=True, help="main's call tables (palm14e_rerun.csv, palm14e.csv)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-cells", type=int, default=0, help="smoke test on the first N cells")
    a = ap.parse_args()
    agg, prv = os.path.join(a.out, "aggregate"), os.path.join(a.out, "private")
    os.makedirs(agg, exist_ok=True)
    os.makedirs(prv, exist_ok=True)
    T = load_targets(a.bed)
    t3 = next(T[T.kind == "indel"].itertuples(index=False))
    tchr = dict(zip(T.tid, T.chr))

    R = pd.read_pickle(a.cache)
    cells = np.sort(R.cell.unique())
    if a.max_cells:
        cells = cells[:a.max_cells]
        R = R[R.cell.isin(set(cells))]
    R = R.drop_duplicates(["cell", "qname", "kind", "pos", "tid"]).reset_index(drop=True)
    E = indel_entries(a.cells, set(cells), t3, a.ref)
    R = R.merge(E, on=["cell", "qname", "pos", "kind", "tid"], how="left")
    for c in ("b20", "b0"):
        R[c] = [e if isinstance(e, tuple) else () for e in R[c]]
    R["chr"] = R.tid.map(tchr)
    R["win"] = (R.pos // WIN) * WIN
    print("records", len(R), "cells", len(cells), "indel-site records", len(E), flush=True)

    def main_grouping(D):
        # awk split in make_consensus.sh: n = all q20 records of (cell, BX, chr, window)
        size = D.drop_duplicates(["cell", "qname", "kind", "pos"]).groupby(["cell", "bx", "chr", "win"]).size()
        return D.join(size.rename("n"), on=["cell", "bx", "chr", "win"])

    P = R[R.kind == "prim"].drop_duplicates(["cell", "qname", "pos", "tid"]).reset_index(drop=True)
    hybrid.SCHEMES = {"T3": (3, True)}
    P = hybrid.assign(P)
    P["n"] = P.groupby(["cell", "tid", "cl_T3"]).cl_T3.transform("size")
    D = {("all", "main"): main_grouping(R), ("prim", "main"): main_grouping(P.drop(columns="n")),
         ("prim", "hybrid"): P}
    keycols = {"main": ["bx", "chr", "win"], "hybrid": ["cl_T3"]}
    grid = pd.MultiIndex.from_product([cells, T.tid], names=["cell", "tid"]).to_frame(index=False)

    tabs, res, ps = {}, {}, {}
    for rung, (recs, grp, nmin, rule, site, caller) in RUNGS.items():
        ck = (recs, grp, nmin, rule, site)
        if ck not in tabs:
            tabs[ck] = cell_table(group_labels(D[(recs, grp)], keycols[grp], nmin, rule, t3, site), grid)
        tab = tabs[ck].copy()
        if caller == "main":
            res[rung] = call_main_exact(tab)
        else:
            res[rung], ps[rung] = call_spec_all(tab, T)
        print(rung, "done", flush=True)

    # ---------------- A0 tie S0 to main's real output on this subset
    tmap = {(r.chr, r.start, r.alt): r.tid for r in T.itertuples()}
    cr = pd.read_csv(a.counts).rename(columns={"barcode": "cell", "MIS": "OTH"})
    cr["tid"] = [tmap[(c, s, x)] for c, s, x in zip(cr.chr, cr.start, cr.alt)]
    cr = grid.merge(cr[["cell", "tid"] + LABELS], on=["cell", "tid"], how="left").fillna(0)
    s0 = res["S0"].set_index(["cell", "tid"])
    cj = cr.set_index(["cell", "tid"])[LABELS].join(s0[LABELS], rsuffix="_emu")
    covered = (cj[LABELS].sum(axis=1) + cj[[c + "_emu" for c in LABELS]].sum(axis=1)) > 0
    a0 = []
    for tid in list(T.tid) + ["all"]:
        s = cj[covered] if tid == "all" else cj[covered & (cj.index.get_level_values("tid") == tid)]
        a0.append({"check": "S0 counts vs counts_rerun.csv", "tid": tid, "rows": len(s),
                   "identical": float((s[LABELS].values == s[[c + "_emu" for c in LABELS]].values).all(axis=1).mean()),
                   **{f"{c}_main": int(s[c].sum()) for c in LABELS}, **{f"{c}_emu": int(s[c + "_emu"].sum()) for c in LABELS}})
    port = call_main_exact(cr.copy()).set_index(["cell", "tid"]).genotype
    for path in a.calls:
        b = pd.read_csv(path)
        b["tid"] = [tmap.get((c, s, x)) for c, s, x in zip(b.chr, b.start, b.alt)]
        b = b.rename(columns={"barcode": "cell"}).dropna(subset=["tid"])
        b = b[b.cell.isin(set(cells))].set_index(["cell", "tid"]).genotype.fillna("NA")
        for name, ours in [("S0 emulation", s0.genotype), ("caller port on counts_rerun", port)]:
            o = ours.reindex(b.index)
            a0.append({"check": f"{name} calls vs {os.path.basename(path)}", "tid": "all", "rows": len(b),
                       "identical": float((o == b).mean()), "agree": int((o == b).sum())})
    pd.DataFrame(a0).to_csv(os.path.join(agg, "A0_s0_tie.csv"), index=False)

    # ---------------- A1-A3 rung comparisons
    a1, a2, a3, flips = [], [], [], []
    for x, y in PAIRS:
        compare(x, y, res[x], res[y], T, a1, a2, a3, flips)
    a1 = pd.DataFrame(a1)
    st = []
    for rung, r in res.items():
        for t in T.itertuples(index=False):
            s = r[r.tid == t.tid]
            c = s[s[LABELS].sum(axis=1) > 0]
            st.append({"rung": rung, "tid": t.tid, "covered_rows": len(c),
                       **{f"umis_{l}": int(s[l].sum()) for l in LABELS},
                       "median_umis_per_covered_row": float(c[LABELS].sum(axis=1).median()) if len(c) else 0,
                       **{f"calls_{g}": int((s.genotype == g).sum()) for g in ("MUT", "WT")},
                       "calls_NA_covered": int((c.genotype == "NA").sum()),
                       "single_umi_MUT_calls": int(((s.genotype == "MUT") & (s.MUT == 1)).sum()),
                       "q_hat": ps[rung][t.tid][0] if rung in ps else np.nan,
                       "P_s": ps[rung][t.tid][1] if rung in ps else np.nan})
    pd.DataFrame(st).to_csv(os.path.join(agg, "A1_rung_state.csv"), index=False)
    a1.to_csv(os.path.join(agg, "A1_rung_flips.csv"), index=False)
    pd.DataFrame(a2).to_csv(os.path.join(agg, "A2_transitions.csv"), index=False)
    pd.DataFrame(a3).to_csv(os.path.join(agg, "A3_hard_flip_kmo.csv"), index=False)
    if flips:
        pd.concat(flips, ignore_index=True).to_csv(os.path.join(prv, "hard_flips.csv"), index=False)

    # ---------------- A4 per cell S0 -> S6
    loc = loci(T)
    j = res["S0"].set_index(["cell", "tid"])[LABELS + ["genotype"]].join(
        res["S6"].set_index(["cell", "tid"])[LABELS + ["genotype"]], rsuffix="_6").reset_index()
    j["changed"] = j.genotype != j.genotype_6
    j["hard"] = j.genotype.isin(HARD) & j.genotype_6.isin(HARD) & j.changed
    j["locus"] = j.tid.map(loc)
    # first rung at which each S0 -> S6 hard flip departs from its S0 call
    chain = pd.concat({r: res[r].set_index(["cell", "tid"]).genotype for r in RUNGS if r not in ("S0+6", "S4n")}, axis=1)
    first = chain.ne(chain.S0, axis=0).idxmax(axis=1).rename("first_rung")
    j = j.join(first, on=["cell", "tid"])
    pc = j.groupby("cell").agg(targets_changed=("changed", "sum"), hard_flips=("hard", "sum"),
                               hard_loci=("locus", lambda s: s[j.loc[s.index, "hard"]].nunique()),
                               umis_S0=("WT", "sum"))
    pc["umis_S0"] += j.groupby("cell")[["MUT", "OTH"]].sum().sum(axis=1)
    faults = pc[pc.hard_loci >= 2]
    jf = j[j.cell.isin(faults.index) & j.hard]
    a4 = [{"metric": "cells", "value": len(pc)},
          {"metric": "cells changing any call", "value": int((pc.targets_changed > 0).sum())},
          {"metric": "cells with a hard flip", "value": int((pc.hard_flips > 0).sum())},
          {"metric": "cells with hard flips at >= 2 independent loci", "value": len(faults)},
          {"metric": "median S0 UMIs per cell (cells with UMIs)", "value": float(pc.umis_S0[pc.umis_S0 > 0].median())},
          {"metric": "median S0 UMIs per cell (fault cells)", "value": float(faults.umis_S0.median()) if len(faults) else np.nan}]
    for n, c in pc.targets_changed.value_counts().sort_index().items():
        a4.append({"metric": f"cells with {n} targets changed", "value": int(c)})
    for (fr, d), c in j[j.hard].assign(d=lambda x: x.genotype + ">" + x.genotype_6).groupby(["first_rung", "d"]).size().items():
        a4.append({"metric": f"S0->S6 hard flips first departing at {fr}, {d}", "value": int(c)})
    for (fr, d), c in jf.assign(d=lambda x: x.genotype + ">" + x.genotype_6).groupby(["first_rung", "d"]).size().items():
        a4.append({"metric": f"fault-cell hard flips first departing at {fr}, {d}", "value": int(c)})
    for pat, c in jf.groupby("cell").tid.agg(lambda s: "+".join(sorted(s))).value_counts().items():
        a4.append({"metric": f"fault cells flipping at {pat}", "value": int(c)})
    pd.DataFrame(a4).to_csv(os.path.join(agg, "A4_per_cell.csv"), index=False)
    pd.DataFrame({"tid": T.tid, "locus": T.tid.map(loc)}).to_csv(os.path.join(agg, "targets_loci.csv"), index=False)
    pc.to_csv(os.path.join(prv, "per_cell_S0_S6.csv"))
    fc = j[j.cell.isin(faults.index)].join(chain.add_prefix("call_"), on=["cell", "tid"])
    fc.to_csv(os.path.join(prv, "fault_cells.csv"), index=False)

    # ---------------- A5 error table
    c0 = res["S0"]
    raw = (c0.groupby("cell").OTH.sum() / (c0.groupby("cell")[LABELS].sum().sum(axis=1) + 1e-8) / 2).rename("raw_err")
    c0 = c0.join(raw, on="cell")
    a5 = []
    for t in T.itertuples(index=False):
        s = c0[(c0.tid == t.tid) & (c0[LABELS].sum(axis=1) > 0)]
        a5.append({"tid": t.tid, "kind": t.kind, "cells_covered": len(s),
                   "main_err_median": s.raw_err.median(), "main_err_q25": s.raw_err.quantile(.25),
                   "main_err_q75": s.raw_err.quantile(.75), "main_err_max": s.raw_err.max(),
                   "cells_err_zero_to_min_fallback": int((s.raw_err == 0).sum()),
                   "cells_at_floor": int((s.err <= EPS).sum()),
                   "main_err_used_median": float(np.maximum(s.err, EPS).median()) if len(s) else np.nan,
                   "v2_q_hat_S6": ps["S6"][t.tid][0], "v2_P_s_S6": ps["S6"][t.tid][1],
                   "v2_q_hat_S0+6": ps["S0+6"][t.tid][0], "v2_P_s_S0+6": ps["S0+6"][t.tid][1]})
    pd.DataFrame(a5).to_csv(os.path.join(agg, "A5_error_table.csv"), index=False)
    for rung, r in res.items():
        r.to_csv(os.path.join(prv, f"calls_{rung}.csv.gz"), index=False)
    print("done", flush=True)


if __name__ == "__main__":
    main()
