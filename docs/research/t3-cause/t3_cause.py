#!/usr/bin/env python3
"""PALM09E T3 cause (#73, map #72): germline, artefact or B-cell over-dispersion.

Implements docs/research/t3-cause-preregistration.md (commit 8b8e330, amended before any data look). Output is aggregate only: targets as
T1..Tn, no barcodes, gene names or coordinates.

Usage:
  t3_cause.py --selftest
  t3_cause.py --sample S --target T3 --v3-dir DIR --qc QC.txt --negatives NEG.csv --metadata META.csv
              --key KEY.txt --bed BED --bam BAM --ref REF.fa --out OUTDIR
"""
import argparse
import csv
import json
import os
import re
from collections import defaultdict

import numpy as np
from scipy import optimize, stats

BARE = re.compile(r"[ACGTN]{16}")
MIN_MAPQ, MIN_MOL_READS, WIN, MIN_BG_READS = 20, 3, 100, 200
COMP = str.maketrans("ACGT", "TGCA")


def bare(c):
    return BARE.search(c).group(0)


def wilson(x, n, z=1.96):
    if n == 0:
        return float("nan"), float("nan")
    p = x / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return float(mid - half), float(mid + half)


# ---------------------------------------------------------------- count level
def bb_test(m, n):
    """Beta-binomial vs binomial LRT over cells (m_i of n_i). Returns r_hat, rho_hat, LR, boundary-corrected p."""
    m, n = np.asarray(m, int), np.asarray(n, int)
    if n.sum() == 0:
        return dict(r=float("nan"), rho=float("nan"), LR=float("nan"), p=float("nan"))
    r0 = m.sum() / n.sum()
    if r0 in (0.0, 1.0):
        return dict(r=float(r0), rho=0.0, LR=0.0, p=1.0)
    l_bin = float(stats.binom.logpmf(m, n, r0).sum())

    def nll(x):
        r, rho = 1 / (1 + np.exp(-x[0])), 1 / (1 + np.exp(-x[1]))
        s = (1 - rho) / rho
        return -float(stats.betabinom.logpmf(m, n, r * s, (1 - r) * s).sum())

    best = None
    for rho0 in (1e-4, 1e-3, 1e-2, 0.05, 0.2, 0.5):
        res = optimize.minimize(nll, [np.log(r0 / (1 - r0)), np.log(rho0 / (1 - rho0))], method="Nelder-Mead",
                                options=dict(xatol=1e-8, fatol=1e-10, maxiter=20000))
        if best is None or res.fun < best.fun:
            best = res
    l_bb = -best.fun
    lr = max(0.0, 2 * (l_bb - l_bin))
    r, rho = 1 / (1 + np.exp(-best.x[0])), 1 / (1 + np.exp(-best.x[1]))
    if l_bb <= l_bin:
        r, rho = r0, 0.0
    return dict(r=float(r), rho=float(rho), LR=float(lr), p=float(0.5 * stats.chi2.sf(lr, 1)) if lr > 0 else 1.0)


def homogeneity(mut, tot, seed=73, draws=100_000):
    """Pearson chi2 homogeneity of MUT rate across groups; Monte Carlo p (fixed margins) if any expected < 5."""
    mut, tot = np.asarray(mut, int), np.asarray(tot, int)
    wt = tot - mut
    table = np.stack([mut, wt])
    exp = table.sum(1, keepdims=True) * table.sum(0, keepdims=True) / table.sum()

    def chi(t):
        with np.errstate(divide="ignore", invalid="ignore"):
            return float(np.nansum((t - exp) ** 2 / exp))

    obs = chi(table)
    if (exp >= 5).all():
        return dict(chi2=obs, p=float(stats.chi2.sf(obs, len(mut) - 1)), method="asymptotic")
    rng = np.random.default_rng(seed)
    sims = rng.multivariate_hypergeometric(tot, int(mut.sum()), size=draws)
    stat = np.array([chi(np.stack([s, tot - s])) for s in sims])
    return dict(chi2=obs, p=float((1 + (stat >= obs - 1e-12).sum()) / (1 + draws)), method="monte_carlo")


def rate_row(name, d):
    x, u = int(d["MUT"].sum()), int(d["n"].sum())
    lo, hi = wilson(x, u)
    cells = len(d)
    cm = int((d["MUT"] >= 1).sum())
    return dict(group=name, cells=cells, umis=u, mut_umis=x, r=x / u if u else float("nan"), r_lo=lo, r_hi=hi,
                cells_with_mut=cm, frac_cells_with_mut=cm / cells if cells else float("nan"),
                cells_mut_ge2=int((d["MUT"] >= 2).sum()))


# ---------------------------------------------------------------- read level
def read_level(bam_path, ref_path, chrom, pos0, ref, alt, group_of):
    import pysam
    fa = pysam.FastaFile(ref_path)
    lo, hi = max(0, pos0 - WIN), pos0 + WIN + 1
    refseq = fa.fetch(chrom, lo, hi).upper()
    assert refseq[pos0 - lo] == ref, f"reference base at T3 is {refseq[pos0 - lo]}, BED REF is {ref}"
    ctx = refseq[pos0 - lo - 1: pos0 - lo + 2]
    t3 = defaultdict(list)          # group -> list of (allele, reverse, dist_to_end, bq)
    mol = defaultdict(lambda: [0, 0])  # (group, CB, UB) -> [reads, alt reads]
    bg = defaultdict(lambda: [0, 0])   # ref position (0-based) -> [ref-or-alt reads, alt reads], Neg reads only
    other_allele = defaultdict(int)
    for rd in pysam.AlignmentFile(bam_path).fetch(chrom, pos0, pos0 + 1):
        if rd.is_unmapped or rd.is_secondary or rd.is_supplementary or rd.mapping_quality < MIN_MAPQ:
            continue
        if not rd.has_tag("CB"):
            continue
        g = group_of.get(bare(rd.get_tag("CB")))
        if g not in ("neg", "leuk"):
            continue
        qs, qq = rd.query_sequence, rd.query_qualities
        pairs = {r: q for q, r in rd.get_aligned_pairs(matches_only=True) if lo <= r < hi}
        q = pairs.get(pos0)
        if q is None:
            other_allele[(g, "del_or_skip")] += 1
            continue
        b = qs[q].upper()
        if b not in (ref, alt):
            other_allele[(g, b)] += 1
            continue
        d = min(q - rd.query_alignment_start, rd.query_alignment_end - 1 - q)
        t3[g].append((b == alt, rd.is_reverse, d, qq[q] if qq is not None else -1))
        if rd.has_tag("UB"):
            m = mol[(g, rd.get_tag("CB"), rd.get_tag("UB"))]
            m[0] += 1
            m[1] += b == alt
        if g == "neg":
            for r, qp in pairs.items():
                if r == pos0 or refseq[r - lo] != ref:
                    continue
                bb = qs[qp].upper()
                if bb in (ref, alt):
                    bg[r][0] += 1
                    bg[r][1] += bb == alt

    out = {"context_trinucleotide_match_count": None}
    for g in ("neg", "leuk"):
        a = np.array(t3[g], dtype=float).reshape(-1, 4)
        isalt, rev, dist, bq = a[:, 0] == 1, a[:, 1] == 1, a[:, 2], a[:, 3]
        n_alt, n_ref = int(isalt.sum()), int((~isalt).sum())
        tab = [[int((isalt & ~rev).sum()), int((isalt & rev).sum())], [int((~isalt & ~rev).sum()), int((~isalt & rev).sum())]]
        fplus_alt = tab[0][0] / n_alt if n_alt else float("nan")
        fplus_ref = tab[1][0] / n_ref if n_ref else float("nan")
        p1 = float(stats.fisher_exact(tab)[1]) if n_alt and n_ref else float("nan")
        if n_alt and n_ref:
            p2 = float(stats.mannwhitneyu(dist[isalt], dist[~isalt], alternative="two-sided").pvalue)
            med_a, med_r = float(np.median(dist[isalt])), float(np.median(dist[~isalt]))
        else:
            p2, med_a, med_r = float("nan"), float("nan"), float("nan")
        out[g] = dict(reads=n_alt + n_ref, alt_reads=n_alt, ref_reads=n_ref,
                      raw_alt_rate=n_alt / (n_alt + n_ref) if n_alt + n_ref else float("nan"),
                      strand_table_altplus_altminus_refplus_refminus=tab, fplus_alt=fplus_alt, fplus_ref=fplus_ref,
                      A1_delta=fplus_alt - fplus_ref, A1_p=p1,
                      dist_median_alt=med_a, dist_median_ref=med_r,
                      A2_ratio=med_a / med_r if med_r else float("nan"), A2_p=p2,
                      bq_median_alt=float(np.median(bq[isalt])) if n_alt else float("nan"),
                      bq_median_ref=float(np.median(bq[~isalt])) if n_ref else float("nan"),
                      other_alleles={k[1]: v for k, v in other_allele.items() if k[0] == g})
    # A3: molecule purity among ALT-majority molecules
    frac = {g: [v[1] / v[0] for k, v in mol.items() if k[0] == g and v[0] >= MIN_MOL_READS] for g in ("neg", "leuk")}
    altm = {g: np.array([f for f in frac[g] if f >= 0.5]) for g in frac}
    for g in frac:
        h = np.histogram(frac[g], bins=[0, 1e-9, 0.2, 0.5, 0.8, 1.0 + 1e-9])[0] if frac[g] else np.zeros(5, int)
        out[g].update(molecules=len(frac[g]), molecules_alt_majority=int(len(altm[g])),
                      molecules_altfrac_hist_0_lt02_lt05_lt08_le1=[int(x) for x in h],
                      alt_majority_purity_median=float(np.median(altm[g])) if len(altm[g]) else float("nan"))
    if len(altm["neg"]) and len(altm["leuk"]):
        out["A3_diff"] = float(np.median(altm["leuk"]) - np.median(altm["neg"]))
        out["A3_p"] = float(stats.mannwhitneyu(altm["neg"], altm["leuk"], alternative="less").pvalue)
    else:
        out["A3_diff"], out["A3_p"] = float("nan"), float("nan")
    # A4: same substitution background in Neg reads
    rates_all, rates_ctx = [], []
    for r, (n, x) in bg.items():
        if n < MIN_BG_READS:
            continue
        rates_all.append(x / n)
        if refseq[r - lo - 1: r - lo + 2] == ctx:
            rates_ctx.append(x / n)
    t3rate = out["neg"]["raw_alt_rate"]
    p95 = float(np.percentile(rates_all, 95)) if rates_all else float("nan")
    p95c = float(np.percentile(rates_ctx, 95)) if rates_ctx else float("nan")
    out.update(A4_bg_positions=len(rates_all), A4_bg_median=float(np.median(rates_all)) if rates_all else float("nan"),
               A4_bg_p95=p95, A4_ratio=t3rate / p95 if p95 else float("inf"),
               A4_ctx_positions=len(rates_ctx), A4_ctx_p95=p95c,
               A4_ctx_ratio=(t3rate / p95c if p95c else float("inf")) if rates_ctx else float("nan"),
               A4_bg_max=float(max(rates_all)) if rates_all else float("nan"))
    out["context_trinucleotide_match_count"] = len(rates_ctx)
    return out


# ---------------------------------------------------------------- decision
def decide(sig):
    G, A, D, U = sig["G"], sig["A"], sig["D"], sig["U"]
    if G:
        return "germline"
    if A and not D:
        return "artefact"
    if U and not A and not D:
        return "uniform signal (germline with allelic imbalance, or molecule-level site artefact)"
    if D and not A:
        return "over-dispersion, " + ("B-specific" if sig["D2"] else "cell-level, not B-specific")
    return "inconclusive"


def germline_signal(rates, nk_cells):
    groups = ["T", "B", "Leuk"] + (["NK"] if nk_cells >= 20 else [])
    het = all(rates[g]["r_lo"] <= 0.7 and rates[g]["r_hi"] >= 0.3 for g in groups)
    hom = all(rates[g]["r_lo"] >= 0.9 for g in groups)
    return bool(het or hom), groups


def selftest():
    rng = np.random.default_rng(1)
    n = rng.integers(1, 6, size=300)
    p_cell = rng.beta(0.04 * 2, 0.96 * 2, size=300)          # rho = 1/3, strongly over-dispersed
    m = rng.binomial(n, p_cell)
    od = bb_test(m, n)
    m0 = rng.binomial(n, 0.04)
    nod = bb_test(m0, n)
    hm = homogeneity([3, 1, 30], [100, 10, 200])
    print(json.dumps(dict(overdispersed=od, binomial=nod, homogeneity_example=hm), indent=1))
    assert od["p"] < 0.05, "self-test: D1 failed to detect a known over-dispersed group"
    assert decide(dict(G=False, A=False, D=True, U=False, D2=True)) == "over-dispersion, B-specific"
    assert decide(dict(G=False, A=True, D=True, U=False, D2=False)) == "inconclusive"
    print("selftest ok")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--selftest", action="store_true")
    for a in ["sample", "target", "v3-dir", "qc", "negatives", "metadata", "key", "bed", "bam", "ref", "out"]:
        p.add_argument(f"--{a}")
    a = p.parse_args()
    if a.selftest:
        return selftest()
    import pandas as pd
    S, TT = a.sample, a.target
    os.makedirs(a.out, exist_ok=True)
    key = [l.strip() for l in open(a.key) if l.strip()]
    gene = key[int(TT[1:]) - 1]

    v3 = pd.read_csv(f"{a.v3_dir}/{S}.csv")
    tg = pd.read_csv(f"{a.v3_dir}/{S}_targets.csv")
    v3["bc"] = v3["cell"].map(bare)
    qc = {bare(l) for l in open(a.qc) if l.strip()}
    neg = {bare(r["cell"]) for r in csv.DictReader(open(a.negatives))} & qc
    meta = pd.read_csv(a.metadata, low_memory=False)
    meta["bc"] = meta["cell16"].map(bare)
    meta = meta[meta["bc"].isin(qc)]
    sets = dict(
        T=set(meta.loc[meta["aggregated_ct"].astype(str).str.startswith("T_"), "bc"]) & neg,
        NK=set(meta.loc[meta["aggregated_ct"] == "NK_cells", "bc"]) & neg,
        B=set(meta.loc[meta["ct"] == "B_cells", "bc"]) & neg,
        Leuk=set(meta.loc[meta["ct"] == "Leukemic", "bc"]) - neg,
        Amb=set(meta.loc[meta["ct"] == "Ambiguous", "bc"]) - neg,
    )
    sets["T/NK"] = sets["T"] | sets["NK"]
    sets["Neg"] = neg
    sets["Other"] = qc - neg - sets["Leuk"] - sets["Amb"]
    dbl = set(meta.loc[meta["scDblFinder.class"].astype(str) == "doublet", "bc"])

    v3 = v3[v3["bc"].isin(qc)].copy()
    v3["n"] = v3["WT"].astype(int) + v3["MUT"].astype(int)
    v3["T"] = v3["gene"].map(lambda g: f"T{key.index(g) + 1}" if g in key else None)
    t3 = v3[(v3["gene"] == gene) & (v3["n"] >= 1)]

    checks = dict(neg_not_T_NK_B=len(neg - sets["T"] - sets["NK"] - sets["B"]),
                  overlap_T_B=len(sets["T"] & sets["B"]), overlap_NK_B=len(sets["NK"] & sets["B"]))
    tn = t3[t3["bc"].isin(neg)]
    c = tn["MUT"].sum() / tn["n"].sum()
    c_run = float(tg.loc[tg["gene"] == gene, "c_s"].iloc[0])
    checks["c_s_recomputed"], checks["c_s_run"] = float(c), c_run
    checks["c_s_match"] = bool(abs(c - c_run) < 1e-12)
    if not checks["c_s_match"]:
        raise SystemExit(f"implementation check failed: {checks}")

    rates = {g: rate_row(g, t3[t3["bc"].isin(sets[g])]) for g in ["T", "NK", "B", "T/NK", "Neg", "Leuk", "Amb", "Other"]}
    D1 = bb_test(tn["MUT"], tn["n"])
    D1_T = bb_test(t3.loc[t3["bc"].isin(sets["T"]), "MUT"], t3.loc[t3["bc"].isin(sets["T"]), "n"])
    D1_L = bb_test(t3.loc[t3["bc"].isin(sets["Leuk"]), "MUT"], t3.loc[t3["bc"].isin(sets["Leuk"]), "n"])
    rb, rt = rates["B"], rates["T/NK"]
    D2_tab = [[rb["cells_with_mut"], rb["cells"] - rb["cells_with_mut"]],
              [rt["cells_with_mut"], rt["cells"] - rt["cells_with_mut"]]]
    D2_p = float(stats.fisher_exact(D2_tab, alternative="greater")[1])
    Uh = homogeneity([rates[g]["mut_umis"] for g in ("T/NK", "B", "Leuk")], [rates[g]["umis"] for g in ("T/NK", "B", "Leuk")])
    G, G_groups = germline_signal(rates, rates["NK"]["cells"])

    # read level
    chrom, start, ref, alt = None, None, None, None
    for l in open(a.bed):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 6 and f[5] == gene:
            chrom, start, ref, alt = f[0], int(f[1]), f[3].upper(), f[4].upper()
    group_of = {b: "neg" for b in neg}
    group_of.update({b: "leuk" for b in sets["Leuk"]})
    rl = read_level(a.bam, a.ref, chrom, start - 1, ref, alt, group_of)
    n_ = rl["neg"]
    A1 = bool(n_["A1_p"] < 0.01 and abs(n_["A1_delta"]) >= 0.25)
    A2 = bool(n_["A2_p"] < 0.01 and n_["A2_ratio"] <= 0.5)
    A3 = bool(rl["A3_p"] < 0.01 and rl["A3_diff"] >= 0.15)
    A4 = bool(rl["A4_ratio"] <= 1)
    D1s, D2s = bool(D1["p"] < 0.05), bool(D2_p < 0.05)
    sig = dict(G=G, A1=A1, A2=A2, A3=A3, A4=A4, A=A1 or A2 or A3 or A4, D1=D1s, D2=D2s, D=D1s or D2s,
               U=bool(Uh["p"] >= 0.05))
    decision = decide(sig)

    # reported, not deciding: cross-target B vs T/NK; MUT-carrying negatives
    cross = []
    for Tn, d in v3[v3["n"] >= 1].groupby("T"):
        b = d[d["bc"].isin(sets["B"])]
        t = d[d["bc"].isin(sets["T/NK"])]
        if len(b) < 20:
            cross.append(dict(target=Tn, B_cells=len(b), note="< 20 covered B cells"))
            continue
        rB, rT = b["MUT"].sum() / b["n"].sum(), (t["MUT"].sum() / t["n"].sum()) if t["n"].sum() else float("nan")
        cross.append(dict(target=Tn, B_cells=len(b), B_umis=int(b["n"].sum()), B_mut=int(b["MUT"].sum()), r_B=float(rB),
                          TNK_cells=len(t), TNK_umis=int(t["n"].sum()), TNK_mut=int(t["MUT"].sum()), r_TNK=float(rT),
                          ratio=float(rB / rT) if rT else float("inf")))
    depth = v3[v3["bc"].isin(neg)].groupby("bc")["n"].sum()
    carriers = set(tn.loc[tn["MUT"] >= 1, "bc"])
    others = set(tn["bc"]) - carriers
    carr = dict(carriers=len(carriers), others=len(others),
                carriers_total_umis_median=float(depth.reindex(list(carriers)).median()) if carriers else float("nan"),
                others_total_umis_median=float(depth.reindex(list(others)).median()) if others else float("nan"),
                carriers_doublet_frac=len(carriers & dbl) / len(carriers) if carriers else float("nan"),
                others_doublet_frac=len(others & dbl) / len(others) if others else float("nan"),
                carriers_by_lineage={g: len(carriers & sets[g]) for g in ("T", "NK", "B")},
                carriers_mut_umis_hist={int(k): int(v) for k, v in tn.loc[tn["MUT"] >= 1, "MUT"].value_counts().items()})
    if carriers and others:
        carr["depth_mw_p"] = float(stats.mannwhitneyu(depth.reindex(list(carriers)), depth.reindex(list(others))).pvalue)

    res = dict(sample=S, target=TT, substitution=f"{ref}>{alt}", checks=checks, rates=rates, G_groups=G_groups,
               D1=D1, D1_T=D1_T, D1_Leuk=D1_L, D2=dict(table_B_TNK_mut_nomut=D2_tab, p=D2_p), U=Uh,
               read_level=rl, signals=sig, decision=decision, cross_target=cross, carriers=carr,
               numpy=np.__version__)
    json.dump(res, open(f"{a.out}/t3_cause_{S}_{TT}.json", "w"), indent=1, default=float)
    pd.DataFrame(rates.values()).to_csv(f"{a.out}/rates_{S}_{TT}.tsv", sep="\t", index=False)
    pd.DataFrame(cross).to_csv(f"{a.out}/cross_target_{S}.tsv", sep="\t", index=False)
    print(json.dumps(dict(signals=sig, decision=decision), indent=1))


if __name__ == "__main__":
    main()
