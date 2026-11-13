#!/usr/bin/env python3
"""Stage 2 of #89: score v3 vs v3-ub against docs/research/confirmation-preregistration.md @ e6f629f, §4–§7, §10.

Reads the stage-1 flags and no-signal list (committed before this stage ran). Contenders share the v3 run's counts:
  v3     θ = min(1, P_s + c)          (c_A held out, c_s in-sample)
  v3-ub  θ = min(1, P_s + CPU(x, n))  (x_A, n_A held out; x_s, n_s in-sample)
Aggregate output only. Usage: score89.py OUTDIR
"""
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2, poisson

import common89 as C
import score_lib as L

CONT = ["v3", "v3-ub"]


def theta(name, P_s, x, n):
    if name == "v3":
        return min(1.0, P_s + (x / n if n else 0.0))
    return min(1.0, P_s + C.cpu(x, n))


def calls_of(rows, th):
    if len(rows) == 0:
        return np.array([], dtype=object)
    _, cl = L.call(rows["MUT"].values, rows["WT"].values, np.full(len(rows), th), 1)
    return cl


def garwood(o, e, conf=0.90):
    a = (1 - conf) / 2
    lo = chi2.ppf(a, 2 * o) / 2 if o > 0 else 0.0
    hi = chi2.ppf(1 - a, 2 * o + 2) / 2
    return (lo / e, hi / e) if e > 0 else (np.nan, np.nan)


def holm(pvals, alpha):
    """Holm step-down; returns dict key -> (adjusted p, rejected)."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m, out, running, stop = len(items), {}, 0.0, False
    for i, (k, p) in enumerate(items):
        adj = min(1.0, (m - i) * p)
        running = max(running, adj)
        rej = (not stop) and p <= alpha / (m - i)
        stop = stop or not rej
        out[k] = (running, rej)
    return out


def med(v):
    return L.seed_summary(v)


def pooled(ts, cont, seed_idx, include):
    """Pooled Poisson-binomial over the included target-seeds of one seed: O, E, p_up, p_low, ratio."""
    o, probs = 0, []
    for r in ts:
        if r["seed"] == seed_idx and include(r):
            o += r[f"O_{cont}"]
            probs.append(r[f"probs_{cont}"])
    pr = np.concatenate(probs) if probs else np.zeros(0)
    e = float(pr.sum())
    return dict(O=o, E=e, p_up=L.poisson_binomial_sf(pr, o), p_low=1.0 - L.poisson_binomial_sf(pr, o + 1),
                ratio=L.ratio(o, e), units=len(probs))


def main():
    out = sys.argv[1]
    C.check_fixed_sets()
    flags = pd.read_csv(f"{out}/stage1_flags.tsv", sep="\t").set_index(["sample", "target"])
    ns = pd.read_csv(f"{out}/stage1_nosignal.tsv", sep="\t").set_index(["sample", "target"])

    ts, tgt_rows, tum_rows, outside = [], [], [], []
    for S in C.PATIENTS:
        d, tg, neg, meta = C.load(S)
        splits = [L.split(neg, s) for s in C.SEEDS]
        for T in tg.index:
            st = C.target_stats(d, tg, T)
            fl = flags.loc[(S, T)]
            assert bool(fl["flag"]) == st["flag"] and abs(fl["c_s"] - st["c_s"]) < 1e-15, f"stage1 drift {S} {T}"
            sub = d[d["T"] == T]
            ng = sub[sub["set"] == "neg"]
            ngc = ng[ng["n"] >= 1]
            cov = sub[sub["n"] >= 1]
            lk, am = cov[cov["set"] == "leuk"], cov[cov["set"] == "amb"]
            row = dict(sample=S, target=T, in_pool=T in C.POOL[S], tumour_evaluable=T in C.TUMOUR.get(S, []),
                       clonal=T in C.CLONAL.get(S, []), flagged=st["flag"], flag_q_clause_only=st["flag_q_clause_only"],
                       P_s=st["P_s"], c_s=st["c_s"], cpu_s=st["cpu_s"], x_s=st["x_s"], n_s=st["n_s"], q_T=st["q_T"],
                       neg_covered=len(ngc), leuk_covered=len(lk), amb_covered=len(am), neg_umis=int(ng["n"].sum()))
            # held-out (pool targets only)
            if T in C.POOL[S]:
                per = {k: [] for k in ["A_cov", "B_cov", "A_umis", "c_A", "cpu_A"]}
                for si, (A, B) in enumerate(splits):
                    a = ngc[ngc["bc"].isin(A)]
                    b = ngc[ngc["bc"].isin(B)]
                    xA, nA = int(a["MUT"].sum()), int(a["n"].sum())
                    r = dict(sample=S, target=T, seed=si, flagged=st["flag"], nA=nA, xA=xA, B_cov=len(b))
                    for k, v in [("A_cov", len(a)), ("B_cov", len(b)), ("A_umis", nA)]:
                        per[k].append(v)
                    if nA > 0:
                        per["c_A"].append(xA / nA)
                        per["cpu_A"].append(C.cpu(xA, nA))
                    for c in CONT:
                        th = theta(c, st["P_s"], xA, nA)
                        cl = calls_of(b, th)
                        r[f"O_{c}"] = int((cl == "MUT").sum())
                        r[f"probs_{c}"] = L.p_mut(b["n"].values, th, 1, th)
                        r[f"E_{c}"] = float(r[f"probs_{c}"].sum())
                    ts.append(r)
                for k, v in per.items():
                    m_, lo, hi = med(v)
                    row.update({f"{k}_med": m_, f"{k}_lo": lo, f"{k}_hi": hi})
                row["zero_info_seeds"] = int(sum(1 for v in per["A_umis"] if v == 0))
            # tumour side, in-sample, all tumour-evaluable targets (gating decided later) + outside-set deployment cost
            rho = min(1.0, st["P_s"] + st["c_s"])
            for c in CONT:
                th = theta(c, st["P_s"], st["x_s"], st["n_s"])
                row[f"theta_in_{c}"] = th
                cl_l, cl_n, cl_a = calls_of(lk, th), calls_of(ngc, th), calls_of(am, th)
                o_t = int((cl_l == "MUT").sum())
                row[f"O_T_{c}"] = o_t
                row[f"O_in_neg_{c}"] = int((cl_n == "MUT").sum())
                lo, hi = L.wilson(row[f"O_in_neg_{c}"], len(ngc))
                row[f"rate_in_neg_{c}"], row[f"rate_in_neg_lo_{c}"], row[f"rate_in_neg_hi_{c}"] = (
                    row[f"O_in_neg_{c}"] / len(ngc) if len(ngc) else np.nan, lo, hi)
                e_t = float(L.p_mut(lk["n"].values, th, 1, rho).sum())
                row[f"E_T_{c}"] = e_t
                row[f"FDR_{c}"] = e_t / o_t if o_t else float("inf")
                for nm, cl in [("neg", cl_n), ("leuk", cl_l), ("amb", cl_a)]:
                    u = int(pd.isna(cl).sum())
                    lo, hi = L.wilson(u, len(cl))
                    row[f"undet_{nm}_{c}"] = u / len(cl) if len(cl) else np.nan
                    row[f"undet_{nm}_lo_{c}"], row[f"undet_{nm}_hi_{c}"] = lo, hi
                if row["tumour_evaluable"] and not st["flag"]:
                    pi, f = float(ns.loc[(S, T), "pi_hat"]), float(ns.loc[(S, T), "f_hat"])
                    row["pi_hat"], row["f_hat"], row["rho0"] = pi, f, rho
                    row[f"efficiency_{c}"] = L.efficiency(o_t, e_t, pi, f, len(lk), rho)
            if (S, T) == (C.PHASED[0], "T3") or (S in (C.PHASED[0],) and T in C.PHASED[1]):
                th = row["theta_in_v3"]
                e_w = float(L.p_mut(lk["n"].values, th, 1, min(1.0, st["P_s"] + st["cpu_s"])).sum())
                row["FDR_v3_weaknull"] = e_w / row["O_T_v3"] if row["O_T_v3"] else float("inf")
            if row["tumour_evaluable"]:
                row["no_signal"] = (None if st["flag"] else str(ns.loc[(S, T), "no_signal"]) == "True")
            tgt_rows.append(row)

    tdf = pd.DataFrame(tgt_rows)

    # ---------------- §5.1 pooled calibration ----------------
    prim = lambda r: (not r["flagged"]) and r["nA"] > 0  # noqa: E731
    frag = lambda r: (not r["flagged"]) and r["nA"] >= C.FRAG_MIN_NA  # noqa: E731
    pooled_rows, pool_sum = [], {}
    for lab, inc in [("primary", prim), ("fragility", frag)]:
        for c in CONT:
            seeds = [pooled(ts, c, si, inc) for si in range(len(C.SEEDS))]
            for si, x in enumerate(seeds):
                pooled_rows.append(dict(test=lab, contender=c, seed=si, **x))
            res = {k: med([x[k] for x in seeds]) for k in ["O", "E", "p_up", "p_low", "ratio", "units"]}
            passes = [x["p_up"] >= 0.05 and x["ratio"] <= 1.5 for x in seeds]
            pool_sum[(lab, c)] = dict(**{f"{k}_{q}": v for k, vv in res.items() for q, v in zip(["med", "lo", "hi"], vv)},
                                      pass_=bool(res["p_up"][0] >= 0.05 and res["ratio"][0] <= 1.5),
                                      seeds_pass=float(np.mean(passes)))
    psum = pd.DataFrame([dict(test=k[0], contender=k[1], **v) for k, v in pool_sum.items()])

    # zero-information target-seeds
    unflag_ts = [r for r in ts if not r["flagged"]]
    zi = pd.DataFrame([dict(sample=r["sample"], zero=r["nA"] == 0, lt11=r["nA"] < C.FRAG_MIN_NA) for r in unflag_ts])
    zi_pat = zi.groupby("sample").agg(target_seeds=("zero", "size"), zero_info=("zero", "sum"), nA_lt11=("lt11", "sum"))
    zi_frac = float(zi["zero"].mean())

    # ---------------- §5.2 per patient + leave-one-patient-out ----------------
    prows = []
    for S in C.PATIENTS:
        for c in CONT:
            rs, lo_, hi_, lopo = [], [], [], []
            for si in range(len(C.SEEDS)):
                sel = [r for r in ts if r["seed"] == si and r["sample"] == S and prim(r)]
                o, e = sum(r[f"O_{c}"] for r in sel), sum(r[f"E_{c}"] for r in sel)
                rs.append(L.ratio(o, e))
                g = garwood(o, e)
                lo_.append(g[0]), hi_.append(g[1])
                lopo.append(pooled(ts, c, si, lambda r, S=S: prim(r) and r["sample"] != S))
            prows.append(dict(sample=S, contender=c, OE_med=med(rs)[0], OE_lo=med(rs)[1], OE_hi=med(rs)[2],
                              garwood90_lo_med=med(lo_)[0], garwood90_hi_med=med(hi_)[0],
                              O_med=med([sum(r[f"O_{c}"] for r in ts if r["seed"] == si and r["sample"] == S and prim(r))
                                         for si in range(len(C.SEEDS))])[0],
                              E_med=med([sum(r[f"E_{c}"] for r in ts if r["seed"] == si and r["sample"] == S and prim(r))
                                         for si in range(len(C.SEEDS))])[0],
                              lopo_OE_med=med([x["ratio"] for x in lopo])[0], lopo_p_med=med([x["p_up"] for x in lopo])[0],
                              lopo_pass=bool(med([x["p_up"] for x in lopo])[0] >= 0.05
                                             and med([x["ratio"] for x in lopo])[0] <= 1.5)))
    pat = pd.DataFrame(prows)

    # ---------------- §5.3 site screen + §10 per-target calibration ----------------
    site_rows = []
    for (S, T), g in pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("probs")} for r in ts]).groupby(
            ["sample", "target"], sort=False):
        flagged = bool(g["flagged"].iloc[0])
        inc = g[g["nA"] > 0]
        row = dict(sample=S, target=T, flagged=flagged, included_seeds=len(inc))
        for c in CONT:
            ps, rr, rate, rlo, rhi, plow = [], [], [], [], [], []
            for _, r in inc.iterrows():
                o, e = int(r[f"O_{c}"]), float(r[f"E_{c}"])
                ps.append(0.0 if (e == 0 and o > 0) else float(poisson.sf(o - 1, e)) if o > 0 else 1.0)
                plow.append(float(poisson.cdf(o, e)) if e > 0 else 1.0)
                rr.append(L.ratio(o, e))
                lo, hi = L.wilson(o, int(r["B_cov"]))
                rate.append(o / r["B_cov"] if r["B_cov"] else np.nan), rlo.append(lo), rhi.append(hi)
            for k, v in [("O_B", inc[f"O_{c}"].tolist()), ("E_B", inc[f"E_{c}"].tolist()), ("OE", rr),
                         ("site_p", ps), ("site_p_low", plow), ("rate_B", rate)]:
                m_, lo, hi = med(v)
                row.update({f"{k}_{c}_med": m_, f"{k}_{c}_lo": lo, f"{k}_{c}_hi": hi})
            row[f"rate_B_wilson_lo_{c}_med"], row[f"rate_B_wilson_hi_{c}_med"] = med(rlo)[0], med(rhi)[0]
        site_rows.append(row)
    site = pd.DataFrame(site_rows)
    site_flags = {}
    for c in CONT:
        tested = site[(~site["flagged"]) & (site["included_seeds"] > 0)]
        h = holm({(r["sample"], r["target"]): r[f"site_p_{c}_med"] for _, r in tested.iterrows()}, C.ALPHA_HOLM)
        site[f"holm_adj_{c}"] = [h.get((r["sample"], r["target"]), (np.nan, None))[0] for _, r in site.iterrows()]
        site[f"site_exception_{c}"] = [h.get((r["sample"], r["target"]), (np.nan, False))[1] for _, r in site.iterrows()]
        site_flags[c] = [f"{s} {t}" for (s, t), (_, rej) in h.items() if rej]
    site["status"] = np.where(site["flagged"], "flagged (excluded)",
                              np.where(site["included_seeds"] == 0, "no negative data", "tested"))

    # ---------------- §6 tumour side ----------------
    tum = tdf[tdf["tumour_evaluable"]].copy()
    tum["status"] = np.where(tum["flagged"], "flagged",
                             np.where(tum["no_signal"].astype(object) == True, "no-signal", "gated"))  # noqa: E712
    gated = tum[tum["status"] == "gated"]
    for c in CONT:
        tum[f"gate_{c}"] = np.where(tum["status"] == "gated", np.where(tum[f"FDR_{c}"] <= 0.2, "pass", "fail"),
                                    tum["status"])

    def locus(S, T):
        return f"{S} {{T3-T5}}" if (S == C.PHASED[0] and T in C.PHASED[1]) else f"{S} {T}"

    tum["locus"] = [locus(s, t) for s, t in zip(tum["sample"], tum["target"])]
    loci = tum.groupby("locus", sort=False)
    locus_status = {k: ("gated" if (g["status"] == "gated").any() else
                        "no-signal" if (g["status"] == "no-signal").all() else
                        "flagged" if (g["status"] == "flagged").all() else "no-signal/flagged")
                    for k, g in loci}
    n_loci = len(locus_status)
    n_ns_loci = sum(1 for v in locus_status.values() if v == "no-signal")
    gated_loci = [k for k, v in locus_status.items() if v == "gated"]
    k_pat = len({k.split()[0] for k in gated_loci})
    tumour_evaluable = not (n_ns_loci >= 3 or k_pat < 3)

    def loss_sum(df, c):
        tot = 0.0
        for k, g in df.groupby("locus", sort=False):
            tot += float(np.median(g[f"O_T_{c}"])) if len(g) > 1 else float(g[f"O_T_{c}"].iloc[0])
        return tot

    gl = tum[tum["status"] == "gated"]
    o3, oub = loss_sum(gl, "v3"), loss_sum(gl, "v3-ub")
    loss = 1 - oub / o3 if o3 else float("nan")
    loss_pat = {S: (1 - loss_sum(g, "v3-ub") / loss_sum(g, "v3")) if loss_sum(g, "v3") else float("nan")
                for S, g in gl.groupby("sample", sort=False)}
    tum["loss_target"] = [1 - r["O_T_v3-ub"] / r["O_T_v3"] if r["O_T_v3"] else np.nan for _, r in tum.iterrows()]
    out_set = tdf[(~tdf.set_index(["sample", "target"]).index.isin(
        list(zip(gl["sample"], gl["target"])))) & (~tdf["flagged"]) & (tdf["leuk_covered"] > 0)]
    deploy = dict(targets=len(out_set), O_T_v3=int(out_set["O_T_v3"].sum()), O_T_ub=int(out_set["O_T_v3-ub"].sum()),
                  silenced_targets=int(((out_set["O_T_v3"] > 0) & (out_set["O_T_v3-ub"] == 0)).sum()))
    weak = tum[tum["locus"] == f"{C.PHASED[0]} {{T3-T5}}"][["target", "status", "FDR_v3", "FDR_v3_weaknull"]]
    wk_gated = weak[weak["status"] == "gated"]
    weak_null_flag = bool(len(wk_gated) and (wk_gated["FDR_v3"] <= 0.2).all() and not (wk_gated["FDR_v3_weaknull"] <= 0.2).all())

    # ---------------- §7 verdict ----------------
    res = {}
    for c in CONT:
        cal = pool_sum[("primary", c)]["pass_"]
        tum_ok = bool((gated[f"FDR_{c}"] <= 0.2).all()) if len(gated) else False
        res[c] = dict(pooled_pass=cal, pooled_pass_fragility=pool_sum[("fragility", c)]["pass_"],
                      tumour_gate_pass=tum_ok, tumour_side_evaluable=tumour_evaluable,
                      tumour_gate_fail=sorted(f"{s} {t}" for s, t, f in zip(gated["sample"], gated["target"],
                                                                             gated[f"FDR_{c}"]) if f > 0.2),
                      passes=bool(cal and tum_ok and tumour_evaluable))
    if res["v3"]["passes"]:
        verdict = "v3 holds"
    elif (not res["v3"]["pooled_pass"]) and res["v3-ub"]["passes"] and loss <= 0.05:
        verdict = "v3-ub replaces v3"
    else:
        verdict = "neither confirmed; v3 stays provisional"
    fragile = any(res[c]["pooled_pass"] != res[c]["pooled_pass_fragility"] for c in CONT)
    under = sorted(f"{S} {T}" for S, ts_ in C.CLONAL.items() for T in ts_ if T not in C.TUMOUR.get(S, []))
    ns_list = sorted(f"{s} {t}" for s, t in zip(tum.loc[tum["status"] == "no-signal", "sample"],
                                                tum.loc[tum["status"] == "no-signal", "target"]))
    scope = (f"Per-target precision confirmed on {len(gated_loci)} covered, signal-bearing loci in {k_pat} patients. "
             f"{len(under)} under-covered and {len(ns_list)} no-signal targets are listed. Thin-null loci: PALM14E {{T3–T5}}."
             + (" Weak-null pass: PALM14E {T3–T5}." if weak_null_flag else ""))
    v = dict(verdict=verdict, contenders=res, loss=loss, loss_gate_pass=bool(loss <= 0.05), loss_by_patient=loss_pat,
             O_T_gated_v3=o3, O_T_gated_ub=oub, deployment_outside_gated=deploy,
             fragile=fragile, site_exceptions=site_flags, weak_null_pass=weak_null_flag,
             zero_info_fraction=zi_frac, zero_info_over_10pct=bool(zi_frac > 0.10),
             loci=n_loci, no_signal_loci=n_ns_loci, gated_loci=gated_loci, patients_with_gated_locus=k_pat,
             locus_status=locus_status, no_signal_targets=ns_list, under_covered=under,
             flagged_pool_targets=sorted(f"{s} {t}" for s, t in zip(site.loc[site["flagged"], "sample"],
                                                                    site.loc[site["flagged"], "target"])),
             no_negative_data=sorted(f"{s} {t}" for s, t in zip(site.loc[site["status"] == "no negative data", "sample"],
                                                                site.loc[site["status"] == "no negative data", "target"])),
             scope_line=scope)

    tdf.to_csv(f"{out}/targets.tsv", sep="\t", index=False)
    pd.DataFrame(pooled_rows).to_csv(f"{out}/pooled_by_seed.tsv", sep="\t", index=False)
    psum.to_csv(f"{out}/pooled.tsv", sep="\t", index=False)
    pat.to_csv(f"{out}/patients.tsv", sep="\t", index=False)
    site.to_csv(f"{out}/site_screen.tsv", sep="\t", index=False)
    tum.to_csv(f"{out}/tumour.tsv", sep="\t", index=False)
    zi_pat.to_csv(f"{out}/zero_info.tsv", sep="\t")
    json.dump(v, open(f"{out}/verdict.json", "w"), indent=1, default=str)
    pd.set_option("display.width", 250)
    print(psum.to_string())
    print(pat.to_string())
    print(json.dumps(v, indent=1, default=str))


if __name__ == "__main__":
    main()
