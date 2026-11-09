#!/usr/bin/env python3
"""Pre-registered scoring (#67, map #60): score `main`, v2, v2 + min2 and v3 against the pass bar of
docs/research/negative-control-preregistration.md (commit 76c3f99), on one sample.

Contenders (pre-registration §2):
  main      main run's counts,  θ_i = max(err_i, baseline_error) per cell,  m_min = 1
  v2        v3 run's counts,    θ = P_s,                                   m_min = 1
  v2+min2   v3 run's counts,    θ = P_s,                                   m_min = 2
  v3        v3 run's counts,    θ = min(1, P_s + c)  (c_A held out, c_s in-sample),  m_min = 1

Output is aggregate only: targets as T1..Tn (key file on the HPC), no barcodes, no gene names.

Usage: score.py --sample S --v3-dir DIR --main-dir DIR --qc QC.txt --negatives NEG.csv --metadata META.csv
                --key KEY.txt --clonal T5,T3,T4 --out OUTDIR [--baseline-error 1e-4]
"""
import argparse
import csv
import json
import os
import re

import numpy as np
import pandas as pd

import score_lib as L

CONTENDERS = [("main", 1), ("v2", 1), ("v2+min2", 2), ("v3", 1)]  # name, m_min
NOISE = 40  # PALM14E run-to-run noise floor, calls
BARE = re.compile(r"[ACGTN]{16}")


def bare(c):
    return BARE.search(c).group(0)


def args_():
    p = argparse.ArgumentParser()
    for a in ["sample", "v3-dir", "main-dir", "qc", "negatives", "metadata", "key", "clonal", "out"]:
        p.add_argument(f"--{a}", required=True)
    p.add_argument("--baseline-error", type=float, default=1e-4)
    return p.parse_args()


def same_calls(a, b):
    return bool((pd.Series(a).fillna("NA").values == pd.Series(b).fillna("NA").values).all())


def calib(rows, theta, m_min):
    """O, E (nominal), p_cal, O/E, observed rate and Wilson CI over covered rows, under rule (θ, m_min)."""
    _, calls = L.call(rows["MUT"].values, rows["WT"].values, theta, m_min)
    o = int((calls == "MUT").sum())
    pm = L.p_mut(rows["n"].values, theta, m_min, theta)
    e = float(pm.sum())
    lo, hi = L.wilson(o, len(rows))
    return dict(O=o, E=e, p_cal=L.poisson_binomial_sf(pm, o), ratio=L.ratio(o, e),
                rate=o / len(rows) if len(rows) else float("nan"), rate_lo=lo, rate_hi=hi, cells=len(rows))


def undet(calls):
    u = int(pd.isna(calls).sum())
    lo, hi = L.wilson(u, len(calls))
    return (u / len(calls) if len(calls) else float("nan")), lo, hi


def one_umi(rows, neg):
    """One-UMI test on covered QC rows: fraction of n = 1 cells reading MUT, negatives vs the other cells."""
    one = rows[rows["n"] == 1]
    out = {}
    for name, sel in [("neg", one["bc"].isin(neg)), ("other", ~one["bc"].isin(neg))]:
        x, d = int((one.loc[sel, "MUT"] == 1).sum()), int(sel.sum())
        lo, hi = L.wilson(x, d)
        out.update({f"one_umi_{name}_cells": d, f"one_umi_mut_frac_{name}": x / d if d else float("nan"),
                    f"one_umi_{name}_lo": lo, f"one_umi_{name}_hi": hi})
    return out


def main():
    a = args_()
    S = a.sample
    os.makedirs(a.out, exist_ok=True)
    key = [l.strip() for l in open(a.key) if l.strip()]
    clonal = set(a.clonal.split(",")) if a.clonal else set()

    v3full = pd.read_csv(f"{a.v3_dir}/{S}.csv")
    tgt = pd.read_csv(f"{a.v3_dir}/{S}_targets.csv")
    mfull = pd.read_csv(f"{a.main_dir}/{S}.csv", index_col=0)
    v3full["bc"] = v3full["cell"].map(bare)
    mfull["bc"] = mfull["cell"].map(bare)
    qc = {bare(l) for l in open(a.qc) if l.strip()}
    neg = {bare(r["cell"]) for r in csv.DictReader(open(a.negatives))}
    meta = pd.read_csv(a.metadata, low_memory=False)
    meta["bc"] = meta["cell16"].map(bare)
    leuk = set(meta.loc[meta["ct"] == "Leukemic", "bc"])
    amb = set(meta.loc[meta["ct"] == "Ambiguous", "bc"])

    # ---- implementation checks (abort on failure) ----
    checks = {}
    _, recalled = L.call(v3full["MUT"].values, v3full["WT"].values, v3full["P_null"].values, 1)
    checks["v3_recall_run_P_null_identical"] = same_calls(recalled, v3full["genotype"])
    mfull["theta"] = np.maximum(mfull["err"].values, a.baseline_error)
    p_imp, cl_imp = L.call(mfull["MUT"].values, mfull["WT"].values, mfull["theta"].values, 1)
    p_run = mfull["p"].values.astype(float)
    # main's nbinom p is NaN when m = 0 (the imported call gives 1 there, k > 0); neither is MUT. Compare on m >= 1.
    m1 = mfull["MUT"].values >= 1
    nan_same = bool(np.isnan(p_run[~m1]).all() and (np.isnan(p_imp[m1]) == np.isnan(p_run[m1])).all())
    fin = m1 & ~np.isnan(p_imp) & ~np.isnan(p_run)
    checks["main_p_nan_pattern_identical"] = nan_same
    checks["main_p_max_abs_diff"] = float(np.max(np.abs(p_imp[fin] - p_run[fin]))) if fin.any() else 0.0
    checks["main_recall_genotype_identical"] = same_calls(cl_imp, mfull["genotype"])  # beyond the pre-registration
    checks["main_rows_unique_cell_gene"] = bool(~mfull.duplicated(["bc", "gene"]).any())
    checks["neg_not_in_qc"] = len(neg - qc)

    v3 = v3full[v3full["bc"].isin(qc)].copy()
    mn = mfull[mfull["bc"].isin(qc)].copy()
    for d in (v3, mn):
        d["n"] = d["WT"].astype(int) + d["MUT"].astype(int)
        d["T"] = d["gene"].map(lambda g: f"T{key.index(g) + 1}" if g in key else None)
        assert d["T"].notna().all(), "target outside the key"
    tid2T = dict(zip(v3["target_id"], v3["T"]))
    max_c_diff, ok = 0.0, True
    for _, r in tgt.iterrows():
        sub = v3[v3["target_id"] == r["target_id"]]
        sn = sub[sub["bc"].isin(neg)]
        c, _, _ = L.pooled_rate(sn["MUT"].values, sn["WT"].values)
        max_c_diff = max(max_c_diff, abs(c - r["c_s"]))
        _, cl = L.call(sub["MUT"].values, sub["WT"].values, np.full(len(sub), min(1.0, r["P_s"] + c)), 1)
        ok &= same_calls(cl, sub["genotype"])
    checks["v3_max_abs_c_s_diff_qc_vs_run"] = max_c_diff
    checks["v3_recall_qc_recomputed_c_s_identical"] = ok
    checks["qc_cells_in_v3"] = int(v3["bc"].nunique())
    checks["qc_cells_in_main"] = int(mn["bc"].nunique())
    checks["numpy"] = np.__version__
    json.dump(checks, open(f"{a.out}/checks_{S}.json", "w"), indent=1)
    if not (checks["v3_recall_run_P_null_identical"] and ok and max_c_diff < 1e-12 and nan_same
            and checks["main_p_max_abs_diff"] < 1e-9 and checks["main_rows_unique_cell_gene"]):
        raise SystemExit(f"implementation check failed: {checks}")

    neg_sorted = sorted(neg & qc)
    splits = [L.split(neg_sorted, s) for s in L.SEEDS]
    v3c, mnc = v3[v3["n"] >= 1], mn[mn["n"] >= 1]

    rows, onum = [], []
    for tid in tgt["target_id"]:
        T = tid2T.get(tid)
        if T is None:
            continue
        P_s = float(tgt.loc[tgt["target_id"] == tid, "P_s"].iloc[0])
        src = {"v3": v3c[v3c["target_id"] == tid], "main": mnc[mnc["T"] == T]}
        part = {k: dict(neg=d[d["bc"].isin(neg)], leuk=d[d["bc"].isin(leuk)], amb=d[d["bc"].isin(amb)])
                for k, d in src.items()}
        t_neg = part["v3"]["neg"]
        c_s, _, _ = L.pooled_rate(t_neg["MUT"].values, t_neg["WT"].values)
        c_main, _, _ = L.pooled_rate(part["main"]["neg"]["MUT"].values, part["main"]["neg"]["WT"].values)
        q_T, _, _ = L.pooled_rate(part["v3"]["leuk"]["MUT"].values, part["v3"]["leuk"]["WT"].values)
        flag = c_s >= 0.25 or (q_T > 0 and c_s >= 0.5 * q_T)
        is_clonal = T in clonal and not flag
        for k in src:
            o = one_umi(src[k], neg)
            onum.append(dict(sample=S, target=T, counts=k, **o))

        # held-out evaluability (v3 counts, all negatives)
        ev, nA_umis = [], []
        for A, B in splits:
            nA = int(t_neg.loc[t_neg["bc"].isin(A), "n"].sum())
            nB = int(t_neg["bc"].isin(B).sum())
            ev.append(nA >= L.MIN_A_UMIS and nB >= L.MIN_B_CELLS)
            nA_umis.append(nA)
        n_ev = int(sum(ev))
        evaluable = n_ev >= L.MIN_EVAL_SEEDS

        # sensitivity mixture fits: v3 counts (shared by v2, v2+min2, v3) and main counts
        fit = {}
        lv = part["v3"]["leuk"]
        rho_v3 = min(1.0, P_s + c_s)
        fit["v3"] = (L.em_mixture(lv["MUT"].values, lv["n"].values, rho_v3) if len(lv) else (np.nan, np.nan),
                     rho_v3)
        lm = part["main"]["leuk"]
        rho_m = np.minimum(1.0, lm["theta"].values + c_main)
        fit["main"] = (L.em_mixture(lm["MUT"].values, lm["n"].values, rho_m) if len(lm) else (np.nan, np.nan),
                       float(np.median(rho_m)) if len(lm) else np.nan)

        def theta_of(name, d, c):
            if name == "main":
                return d["theta"].values
            return np.full(len(d), min(1.0, P_s + (c if name == "v3" else 0.0)))

        base = {}
        for name, m_min in CONTENDERS:
            cs = "main" if name == "main" else "v3"
            P = part[cs]
            # in-sample
            ins = calib(P["neg"], theta_of(name, P["neg"], c_s), m_min)
            th_l = theta_of(name, P["leuk"], c_s)
            _, cl_leuk = L.call(P["leuk"]["MUT"].values, P["leuk"]["WT"].values, th_l, m_min)
            _, cl_neg = L.call(P["neg"]["MUT"].values, P["neg"]["WT"].values, theta_of(name, P["neg"], c_s), m_min)
            _, cl_amb = L.call(P["amb"]["MUT"].values, P["amb"]["WT"].values, theta_of(name, P["amb"], c_s), m_min)
            rho_true = (np.minimum(1.0, P["leuk"]["theta"].values + c_main) if name == "main"
                        else np.full(len(P["leuk"]), rho_v3))
            o_t = int((cl_leuk == "MUT").sum())
            e_t = float(L.p_mut(P["leuk"]["n"].values, th_l, m_min, rho_true).sum())
            e_t_nom = float(L.p_mut(P["leuk"]["n"].values, th_l, m_min, th_l).sum())
            fdr = e_t / o_t if o_t else float("inf")
            (pi, f), rho0 = fit[cs]
            eff = L.efficiency(o_t, e_t, pi, f, len(P["leuk"]), rho0) if len(P["leuk"]) else float("nan")
            un = {k: undet(v) for k, v in [("neg", cl_neg), ("leuk", cl_leuk), ("amb", cl_amb)]}
            # held out
            ho = {k: [] for k in ["c_A", "O_B", "E_B", "p_cal", "ratio", "rate", "B_cells", "pass"]}
            for (A, B), e_ok in zip(splits, ev):
                if not e_ok:
                    continue
                c_A, _, _ = L.pooled_rate(t_neg.loc[t_neg["bc"].isin(A), "MUT"].values,
                                          t_neg.loc[t_neg["bc"].isin(A), "WT"].values)
                tb = P["neg"][P["neg"]["bc"].isin(B)]
                r = calib(tb, theta_of(name, tb, c_A), m_min)
                ho["c_A"].append(c_A)
                ho["O_B"].append(r["O"])
                ho["E_B"].append(r["E"])
                for k in ["p_cal", "ratio", "rate"]:
                    ho[k].append(r[k])
                ho["B_cells"].append(r["cells"])
                ho["pass"].append(r["p_cal"] >= 0.05 and r["ratio"] <= 1.5)
            med = {k: L.seed_summary(v) for k, v in ho.items() if k != "pass"}
            c1 = ("pass" if med["p_cal"][0] >= 0.05 and med["ratio"][0] <= 1.5 else "fail") if evaluable \
                else "not evaluable"
            c2 = ("pass" if fdr <= 0.2 else "fail") if T in clonal else "-"
            if flag:  # still scored and reported, excluded from the verdict (pre-registration §4)
                c1, c2 = f"flagged ({c1})", f"flagged ({c2})"
            base[name] = o_t
            rows.append(dict(
                sample=S, target=T, clonal=T in clonal, flagged=flag, contender=name, counts=cs, m_min=m_min,
                P_s=P_s, c_s=c_s, c_s_main=c_main, q_T=q_T,
                theta_median=float(np.median(theta_of(name, P["neg"], c_s))) if len(P["neg"]) else np.nan,
                neg_cells=len(P["neg"]), leuk_cells=len(P["leuk"]), amb_cells=len(P["amb"]),
                eval_seeds=n_ev, evaluable=evaluable, A_umis_med=float(np.median(nA_umis)),
                B_cells_med=med["B_cells"][0],
                c_A_med=med["c_A"][0], c_A_lo=med["c_A"][1], c_A_hi=med["c_A"][2],
                O_B_med=med["O_B"][0], O_B_lo=med["O_B"][1], O_B_hi=med["O_B"][2],
                E_B_med=med["E_B"][0], E_B_lo=med["E_B"][1], E_B_hi=med["E_B"][2],
                ratio_med=med["ratio"][0], ratio_lo=med["ratio"][1], ratio_hi=med["ratio"][2],
                p_cal_med=med["p_cal"][0], p_cal_lo=med["p_cal"][1], p_cal_hi=med["p_cal"][2],
                rate_B_med=med["rate"][0], rate_B_lo=med["rate"][1], rate_B_hi=med["rate"][2],
                seeds_pass=(sum(ho["pass"]) / len(ho["pass"])) if ho["pass"] else float("nan"),
                O_in=ins["O"], E_in=ins["E"], ratio_in=ins["ratio"], p_cal_in=ins["p_cal"],
                rate_in=ins["rate"], rate_in_lo=ins["rate_lo"], rate_in_hi=ins["rate_hi"],
                O_T=o_t, E_T=e_t, E_T_nominal=e_t_nom, FDR=fdr, pi_hat=pi, f_hat=f, rho0=rho0, efficiency=eff,
                eff_pass=(eff >= 0.7) if not np.isnan(eff) else None,
                undet_neg=un["neg"][0], undet_neg_lo=un["neg"][1], undet_neg_hi=un["neg"][2],
                undet_leuk=un["leuk"][0], undet_leuk_lo=un["leuk"][1], undet_leuk_hi=un["leuk"][2],
                undet_amb=un["amb"][0], undet_amb_lo=un["amb"][1], undet_amb_hi=un["amb"][2],
                crit1=c1, crit2=c2,
            ))
        for r in rows:
            if r["sample"] == S and r["target"] == T:
                d = r["O_T"] - base["v2+min2"]
                r["dO_T_vs_v2min2"] = d
                r["dO_T_vs_v3"] = r["O_T"] - base["v3"]
                c3 = ("pass" if d >= NOISE else "tie") if (r["contender"] == "v3" and T in clonal) else "-"
                r["crit3a"] = f"flagged ({c3})" if flag else c3

    out = pd.DataFrame(rows)
    out.to_csv(f"{a.out}/score_{S}.tsv", sep="\t", index=False)
    pd.DataFrame(onum).to_csv(f"{a.out}/one_umi_{S}.tsv", sep="\t", index=False)

    verdict = {}
    for name, _ in CONTENDERS:
        x = out[(out["contender"] == name) & ~out["flagged"]]
        fail1 = sorted(x.loc[x["evaluable"] & (x["crit1"] != "pass"), "target"])
        fail2 = sorted(x.loc[x["clonal"] & (x["crit2"] != "pass"), "target"])
        fail3 = sorted(x.loc[x["clonal"] & (x["crit3a"] != "pass"), "target"]) if name == "v3" else []
        verdict[name] = dict(crit1_fail=fail1, crit2_fail=fail2, crit3a_not_met=fail3,
                             sample_pass=not (fail1 or fail2 or fail3))
    json.dump(verdict, open(f"{a.out}/verdict_{S}.json", "w"), indent=1)
    print(json.dumps(dict(checks=checks, verdict=verdict), indent=1))


if __name__ == "__main__":
    main()
