#!/usr/bin/env python3
"""Exploratory caller-only sweep (#68, map #60): rescore the genotype-v3 run's consensus counts with caller
variants outside the pre-registered verdict, using the pre-registration's computations (score_lib).

Contenders (all from the v3 run's counts; P_s from the run, not re-estimated):
  v2        θ = P_s,                 m_min = 1   (reference)
  v2+min2   θ = P_s,                 m_min = 2   (reference; the 3a baseline)
  v3        θ = P_s + c,             m_min = 1   (reference; c from all T/B/NK negatives)
  v3-ub     θ = P_s + CP95(c),       m_min = 1   (one-sided 95 % Clopper-Pearson upper bound of c)
  v3-tb     θ = P_s + c[T/B only],   m_min = 1   (NK cells dropped from the c estimate; scoring cells unchanged)
  v3-min2   θ = P_s + c,             m_min = 2

Output is aggregate only: targets as T1..Tn (key file on the HPC), no barcodes, no gene names.

Usage: sweep.py --sample S --v3-dir DIR --qc QC.txt --negatives NEG.csv --metadata META.csv --key KEY.txt
                --clonal T5,T3,T7 --out OUTDIR
"""
import argparse
import csv
import json
import os
import re

import numpy as np
import pandas as pd

import score_lib as L

CONTENDERS = [  # name, c source, c bound, m_min
    ("v2", None, "point", 1),
    ("v2+min2", None, "point", 2),
    ("v3", "all", "point", 1),
    ("v3-ub", "all", "cp95", 1),
    ("v3-tb", "tb", "point", 1),
    ("v3-min2", "all", "point", 2),
]
BARE = re.compile(r"[ACGTN]{16}")


def bare(c):
    return BARE.search(c).group(0)


def args_():
    p = argparse.ArgumentParser()
    for a in ["sample", "v3-dir", "qc", "negatives", "metadata", "key", "clonal", "out"]:
        p.add_argument(f"--{a}", required=True)
    return p.parse_args()


def c_of(sub, source, bound, nk):
    """Contamination estimate of a contender from a set of negative rows `sub` at one target."""
    if source is None:
        return 0.0
    if source == "tb":
        sub = sub[~sub["bc"].isin(nk)]
    c, x, d = L.pooled_rate(sub["MUT"].values, sub["WT"].values)
    return L.cp_upper(x, d) if bound == "cp95" else c


def calib(rows, theta, m_min, calls):
    """O, E (nominal), p_cal, O/E, observed rate and Wilson CI over covered rows."""
    o = int((calls == "MUT").sum())
    pm = L.p_mut(rows["n"].values, theta, m_min, theta)
    e = float(pm.sum())
    lo, hi = L.wilson(o, len(rows))
    return dict(O=o, E=e, p_cal=L.poisson_binomial_sf(pm, o), ratio=L.ratio(o, e),
                rate=o / len(rows) if len(rows) else float("nan"), rate_lo=lo, rate_hi=hi, cells=len(rows))


def main():
    a = args_()
    S = a.sample
    os.makedirs(a.out, exist_ok=True)
    key = [l.strip() for l in open(a.key) if l.strip()]
    clonal = set(a.clonal.split(",")) if a.clonal else set()

    full = pd.read_csv(f"{a.v3_dir}/{S}.csv")
    tgt = pd.read_csv(f"{a.v3_dir}/{S}_targets.csv")
    full["bc"] = full["cell"].map(bare)
    qc = {bare(l) for l in open(a.qc) if l.strip()}
    neg = {bare(r["cell"]) for r in csv.DictReader(open(a.negatives))}
    meta = pd.read_csv(a.metadata, low_memory=False)
    meta["bc"] = meta["cell16"].map(bare)
    nk = set(meta.loc[meta["aggregated_ct"] == "NK_cells", "bc"])
    tcell = set(meta.loc[meta["aggregated_ct"].astype(str).str.startswith("T_"), "bc"])
    bcell = set(meta.loc[meta["ct"] == "B_cells", "bc"])
    leuk = set(meta.loc[meta["ct"] == "Leukemic", "bc"])
    amb = set(meta.loc[meta["ct"] == "Ambiguous", "bc"])

    # ---- implementation checks (abort on failure) ----
    checks = {}
    _, recalled = L.call(full["MUT"].values, full["WT"].values, full["P_null"].values, 1)
    same = (pd.Series(recalled).fillna("NA").values == full["genotype"].fillna("NA").values)
    checks["recall_run_P_null_identical"] = bool(same.all())
    checks["neg_not_in_qc"] = len(neg - qc)
    checks["nk_in_neg"] = len(nk & neg & qc)
    checks["nk_not_in_neg"] = len((nk & qc) - neg)
    checks["neg_not_T_B_NK"] = len((neg & qc) - nk - tcell - bcell)
    df = full[full["bc"].isin(qc)].copy()
    df["n"] = df["WT"].astype(int) + df["MUT"].astype(int)
    df["T"] = df["gene"].map(lambda g: f"T{key.index(g) + 1}" if g in key else None)
    assert df["T"].notna().all(), "target outside the key"
    tid2T = dict(zip(df["target_id"], df["T"]))
    max_c_diff = 0.0
    for _, r in tgt.iterrows():
        sub = df[(df["target_id"] == r["target_id"]) & df["bc"].isin(neg)]
        c, _, _ = L.pooled_rate(sub["MUT"].values, sub["WT"].values)
        max_c_diff = max(max_c_diff, abs(c - r["c_s"]))
    checks["max_abs_c_s_diff_qc_vs_run"] = max_c_diff
    # re-call QC rows with the recomputed in-sample P_null
    ok = True
    for _, r in tgt.iterrows():
        sub = df[df["target_id"] == r["target_id"]]
        c, _, _ = L.pooled_rate(sub.loc[sub["bc"].isin(neg), "MUT"].values, sub.loc[sub["bc"].isin(neg), "WT"].values)
        _, cl = L.call(sub["MUT"].values, sub["WT"].values, np.full(len(sub), min(1.0, r["P_s"] + c)), 1)
        ok &= bool((pd.Series(cl).fillna("NA").values == sub["genotype"].fillna("NA").values).all())
    checks["recall_qc_recomputed_c_s_identical"] = ok
    checks["numpy"] = np.__version__
    json.dump(checks, open(f"{a.out}/checks_{S}.json", "w"), indent=1)
    if not (checks["recall_run_P_null_identical"] and ok and max_c_diff < 1e-12):
        raise SystemExit(f"implementation check failed: {checks}")

    neg_sorted = sorted(neg & qc)
    splits = [L.split(neg_sorted, s) for s in L.SEEDS]
    cov = df[df["n"] >= 1]

    rows, groups = [], []
    for tid in tgt["target_id"]:
        T = tid2T.get(tid)
        if T is None:
            continue
        P_s = float(tgt.loc[tgt["target_id"] == tid, "P_s"].iloc[0])
        t_all = cov[cov["target_id"] == tid]
        t_neg = t_all[t_all["bc"].isin(neg)]
        t_leuk = t_all[t_all["bc"].isin(leuk)]
        t_amb = t_all[t_all["bc"].isin(amb)]
        c_s, _, _ = L.pooled_rate(t_neg["MUT"].values, t_neg["WT"].values)
        q_T, _, _ = L.pooled_rate(t_leuk["MUT"].values, t_leuk["WT"].values)
        flag = c_s >= 0.25 or (q_T > 0 and c_s >= 0.5 * q_T)
        rho_true = min(1.0, P_s + c_s)
        pi, f = L.em_mixture(t_leuk["MUT"].values, t_leuk["n"].values, rho_true) if len(t_leuk) else (np.nan, np.nan)

        # contamination by negative subgroup (descriptive)
        for gname, gset in [("T", tcell), ("B", bcell), ("NK", nk), ("T/B", (tcell | bcell)), ("T/B/NK", neg)]:
            g = t_neg[t_neg["bc"].isin(gset)]
            c, x, d = L.pooled_rate(g["MUT"].values, g["WT"].values)
            lo, hi = L.wilson(x, d)
            groups.append(dict(sample=S, target=T, group=gname, cells=len(g), umis=d, mut=x, c=c, lo=lo, hi=hi))

        # held-out evaluability (v3 counts, all negatives)
        ev = []
        for A, B in splits:
            nA = int(t_neg.loc[t_neg["bc"].isin(A), "n"].sum())
            nB = int(t_neg["bc"].isin(B).sum())
            ev.append(nA >= L.MIN_A_UMIS and nB >= L.MIN_B_CELLS)
        n_ev = int(sum(ev))
        evaluable = n_ev >= L.MIN_EVAL_SEEDS

        base = {}
        for name, src, bound, m_min in CONTENDERS:
            c_in = c_of(t_neg, src, bound, nk)
            th = min(1.0, P_s + c_in)
            # in-sample
            _, cl_neg = L.call(t_neg["MUT"].values, t_neg["WT"].values, np.full(len(t_neg), th), m_min)
            ins = calib(t_neg, th, m_min, cl_neg)
            _, cl_leuk = L.call(t_leuk["MUT"].values, t_leuk["WT"].values, np.full(len(t_leuk), th), m_min)
            _, cl_amb = L.call(t_amb["MUT"].values, t_amb["WT"].values, np.full(len(t_amb), th), m_min)
            o_t = int((cl_leuk == "MUT").sum())
            e_t = float(L.p_mut(t_leuk["n"].values, th, m_min, rho_true).sum())
            e_t_nom = float(L.p_mut(t_leuk["n"].values, th, m_min, th).sum())
            fdr = e_t / o_t if o_t else float("inf")
            eff = L.efficiency(o_t, e_t, pi, f, len(t_leuk), rho_true) if len(t_leuk) else float("nan")
            und = {}
            for cname, cl in [("neg", cl_neg), ("leuk", cl_leuk), ("amb", cl_amb)]:
                u = int(pd.isna(cl).sum())
                lo, hi = L.wilson(u, len(cl))
                und[cname] = (u / len(cl) if len(cl) else float("nan"), lo, hi)
            # held out
            ho = {k: [] for k in ["c_A", "O_B", "E_B", "p_cal", "ratio", "rate", "pass"]}
            for (A, B), e_ok in zip(splits, ev):
                if not e_ok:
                    continue
                c_A = c_of(t_neg[t_neg["bc"].isin(A)], src, bound, nk)
                thA = min(1.0, P_s + c_A)
                tb = t_neg[t_neg["bc"].isin(B)]
                _, clB = L.call(tb["MUT"].values, tb["WT"].values, np.full(len(tb), thA), m_min)
                r = calib(tb, thA, m_min, clB)
                ho["c_A"].append(c_A)
                for k in ["O_B", "E_B", "p_cal", "ratio", "rate"]:
                    ho[k].append(r[k.replace("_B", "")] if k in ("O_B", "E_B") else r[k])
                ho["pass"].append(r["p_cal"] >= 0.05 and r["ratio"] <= 1.5)
            med = {k: L.seed_summary(v) for k, v in ho.items() if k != "pass"}
            if not evaluable:
                c1 = "not evaluable"
            else:
                c1 = "pass" if med["p_cal"][0] >= 0.05 and med["ratio"][0] <= 1.5 else "fail"
            c2 = ("pass" if fdr <= 0.2 else "fail") if T in clonal else "-"
            base[name] = o_t
            rows.append(dict(
                sample=S, target=T, clonal=T in clonal, flagged=flag, contender=name, m_min=m_min,
                P_s=P_s, c_s=c_s, c_used=c_in, theta=th, q_T=q_T,
                neg_cells=len(t_neg), leuk_cells=len(t_leuk), amb_cells=len(t_amb),
                eval_seeds=n_ev, evaluable=evaluable,
                c_A_med=med["c_A"][0], c_A_lo=med["c_A"][1], c_A_hi=med["c_A"][2],
                O_B_med=med["O_B"][0], O_B_lo=med["O_B"][1], O_B_hi=med["O_B"][2],
                E_B_med=med["E_B"][0], ratio_med=med["ratio"][0], ratio_lo=med["ratio"][1], ratio_hi=med["ratio"][2],
                p_cal_med=med["p_cal"][0], p_cal_lo=med["p_cal"][1], p_cal_hi=med["p_cal"][2],
                rate_B_med=med["rate"][0], seeds_pass=(sum(ho["pass"]) / len(ho["pass"])) if ho["pass"] else float("nan"),
                O_in=ins["O"], E_in=ins["E"], ratio_in=ins["ratio"], p_cal_in=ins["p_cal"],
                rate_in=ins["rate"], rate_in_lo=ins["rate_lo"], rate_in_hi=ins["rate_hi"],
                O_T=o_t, E_T=e_t, E_T_nominal=e_t_nom, FDR=fdr, pi_hat=pi, f_hat=f, efficiency=eff,
                undet_neg=und["neg"][0], undet_neg_lo=und["neg"][1], undet_neg_hi=und["neg"][2],
                undet_leuk=und["leuk"][0], undet_leuk_lo=und["leuk"][1], undet_leuk_hi=und["leuk"][2],
                undet_amb=und["amb"][0], crit1=c1 if not flag else f"flagged ({c1})", crit2=c2 if not flag else "flagged",
            ))
        for r in rows:
            if r["sample"] == S and r["target"] == T:
                d = r["O_T"] - base["v2+min2"]
                r["dO_T_vs_v2min2"] = d
                r["crit3a"] = ("pass" if d >= 40 else "tie") if (r["contender"].startswith("v3") and T in clonal and not flag) else "-"

    pd.DataFrame(rows).to_csv(f"{a.out}/sweep_{S}.tsv", sep="\t", index=False)
    pd.DataFrame(groups).to_csv(f"{a.out}/c_by_group_{S}.tsv", sep="\t", index=False)

    # sample verdict per contender (exploratory: same rule as the pre-registration, provisional clonal set)
    out = pd.DataFrame(rows)
    verdict = {}
    for name, *_ in CONTENDERS:
        x = out[(out["contender"] == name) & ~out["flagged"]]
        ok1 = (x.loc[x["evaluable"], "crit1"] == "pass").all()
        ok2 = (x.loc[x["clonal"], "crit2"] == "pass").all()
        ok3 = (x.loc[x["clonal"], "crit3a"] == "pass").all() if name.startswith("v3") else True
        verdict[name] = dict(crit1=bool(ok1), crit2=bool(ok2), crit3a=bool(ok3), sample_pass=bool(ok1 and ok2 and ok3))
    json.dump(verdict, open(f"{a.out}/verdict_{S}.json", "w"), indent=1)
    print(json.dumps(dict(checks=checks, verdict=verdict), indent=1))


if __name__ == "__main__":
    main()
