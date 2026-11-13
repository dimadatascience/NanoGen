#!/usr/bin/env python3
"""Ticket #91 (map #72): secondary arm, position-tolerant indel matching, scored against
docs/research/confirmation-preregistration.md @ e6f629f §8 (and §5.1, §5.3, §6.3 by reference).

The primary verdict (#90, c637110) is rule 4: v3 stays provisional, so the §8 call rule is **v3**:
  θ = min(1, P_s + c)   (c_A held out, c_s in-sample), m_min = 1, MUT iff m ≥ 1 and p < 0.1.
Contender = the indel-arm run (prototype c77e147, --indel_window 12), with P_s, c_A, c_s and calls from its own
counts. Baseline = the same rule on the v3 run. Outlier flags and no-signal status are the frozen #89 stage-1
tables (v3 run, §5.4, §6.2). Library: common89 + score_lib @ 23e407c, `call` from genotype-v3@61d8809.
Aggregate output only; targets as T1..Tn. Usage: score91.py Q89_OUT OUTDIR
"""
import json
import math
import sys

import numpy as np
import pandas as pd
from scipy.stats import poisson

import common89 as C
import score_lib as L
from score89 import holm, pooled, calls_of, med

ARMS = ["v3", "proto"]
# §8: the pool's indel targets (10) and the gain set (tumour-evaluable clonal indel targets, before flag/no-signal).
INDEL_POOL = {"PALM08_E": ["T3", "T6", "T7"], "PALM11_E": ["T1"], "PALM14_E": ["T3"],
              "PALM21_E": ["T1", "T2", "T7"], "AML5": ["T2"], "sAML1": ["T11"]}
GAIN_SET = [("PALM08_E", "T3"), ("PALM11_E", "T1"), ("PALM14_E", "T3"), ("sAML1", "T11")]
# §8 tie band: max(ceil(40 * rows_s / 34,979), measured_s); measured_s is computed below from SNV rows.
FLOOR_N, FLOOR_ROWS = 40, 34979


def theta_v3(P_s, x, n):
    return min(1.0, P_s + (x / n if n else 0.0))


def load_arm(S, arm):
    """C.load for either run directory (the v3 path is hard-coded there)."""
    run = "v3" if arm == "v3" else "proto"
    inputs = pd.read_csv(f"{C.W}/q79_full/sheets/inputs.tsv", sep="\t").set_index("sample")
    d = pd.read_csv(f"{C.RUNS}/{run}_{S}/results/{S}/{S}.csv")
    tg = pd.read_csv(f"{C.RUNS}/{run}_{S}/results/{S}/{S}_targets.csv")
    keys = C.bed_keys(inputs.loc[S, "bed"])
    t_of = {}
    for i, (c, s, e, r, a) in enumerate(keys):
        cand = {f"{c}:{p}:{r}>{a}" for p in (s, s + 1, e)}
        hit = [t for t in tg["target_id"] if t in cand]
        assert len(hit) == 1, f"{S}: BED row {i + 1} matches {len(hit)} target ids"
        t_of[hit[0]] = f"T{i + 1}"
    qc = C.read_list(inputs.loc[S, "qc"])
    neg = C.read_list(inputs.loc[S, "negatives"]) & qc
    lab = C.read_labels(C.label_path(S))
    leuk = {b for b, c in lab.items() if c == "Leukemic"} & qc
    amb = {b for b, c in lab.items() if c == "Ambiguous"} & qc
    d["bc"] = d["cell"].map(C.bare)
    d = d[d["bc"].isin(qc)].copy()
    d["n"] = d["WT"].astype(int) + d["MUT"].astype(int)
    d["T"] = d["target_id"].map(t_of)
    d["set"] = np.where(d["bc"].isin(neg), "neg", np.where(d["bc"].isin(leuk), "leuk",
                                                           np.where(d["bc"].isin(amb), "amb", "other")))
    tg["T"] = tg["target_id"].map(t_of)
    return d, tg.set_index("T"), sorted(neg)


def main():
    q89, out = sys.argv[1], sys.argv[2]
    C.check_fixed_sets()
    flags = pd.read_csv(f"{q89}/stage1_flags.tsv", sep="\t").set_index(["sample", "target"])
    ns = pd.read_csv(f"{q89}/stage1_nosignal.tsv", sep="\t").set_index(["sample", "target"])
    pooled89 = pd.read_csv(f"{q89}/pooled.tsv", sep="\t")

    checks, bands, ts, trows = {}, {}, [], []
    for S in C.PATIENTS:
        dv, tgv, neg = load_arm(S, "v3")
        dp, tgp, negp = load_arm(S, "proto")
        assert neg == negp
        key = ["bc", "T"]
        m = dv[key + ["WT", "MUT", "genotype", "set", "n"]].merge(
            dp[key + ["WT", "MUT", "genotype", "n"]], on=key, suffixes=("_v", "_p"), how="outer", indicator=True)
        assert (m["_merge"] == "both").all(), f"{S}: row sets differ"
        indel_T = set(tgv.index[tgv["type"] == "INDEL"])
        assert indel_T == set(tgp.index[tgp["type"] == "INDEL"])
        assert indel_T & set(C.POOL[S]) == set(INDEL_POOL.get(S, [])), f"{S}: indel pool {indel_T}"
        snv = m[~m["T"].isin(indel_T)]
        gv, gp = snv["genotype_v"].fillna("NA"), snv["genotype_p"].fillna("NA")
        snv_diff = int(((snv["WT_v"] != snv["WT_p"]) | (snv["MUT_v"] != snv["MUT_p"]) | (gv != gp)).sum())
        rows_s = len(dv)
        measured = math.ceil(snv_diff * rows_s / len(snv)) if len(snv) else 0
        floor = math.ceil(FLOOR_N * rows_s / FLOOR_ROWS)
        bands[S] = dict(rows_s=rows_s, snv_rows=len(snv), snv_diff=snv_diff, measured_s=measured, floor=floor,
                        band=max(floor, measured))
        # implementation check on the prototype: re-call with its own in-sample P_null reproduces its genotype column
        recall_diff = 0
        for T in tgp.index:
            sub = dp[dp["T"] == T]
            _, cl = L.call(sub["MUT"].values, sub["WT"].values, np.full(len(sub), float(tgp.loc[T, "P_null"])), 1)
            recall_diff += int((pd.Series(cl, index=sub.index).fillna("NA") != sub["genotype"].fillna("NA")).sum())
        checks[S] = dict(recall_diff_proto=recall_diff, indel_targets=sorted(indel_T, key=lambda t: int(t[1:])),
                         indel_rows_call_diff=int(((m["genotype_v"].fillna("NA") != m["genotype_p"].fillna("NA"))
                                                   & m["T"].isin(indel_T)).sum()))

        splits = [L.split(neg, s) for s in C.SEEDS]
        for T in tgv.index:
            fl = flags.loc[(S, T)]
            base = dict(sample=S, target=T, indel=T in indel_T, in_pool=T in C.POOL[S],
                        tumour_evaluable=T in C.TUMOUR.get(S, []), flagged=bool(fl["flag"]))
            if base["tumour_evaluable"] and not base["flagged"]:
                base["no_signal"] = str(ns.loc[(S, T), "no_signal"]) == "True"
            for arm, d, tg in [("v3", dv, tgv), ("proto", dp, tgp)]:
                st = C.target_stats(d, tg, T)
                if arm == "v3":
                    assert abs(st["c_s"] - fl["c_s"]) < 1e-15, f"stage1 drift {S} {T}"
                sub = d[d["T"] == T]
                ngc = sub[(sub["set"] == "neg") & (sub["n"] >= 1)]
                lk = sub[(sub["set"] == "leuk") & (sub["n"] >= 1)]
                th = theta_v3(st["P_s"], st["x_s"], st["n_s"])
                rho = min(1.0, st["P_s"] + st["c_s"])
                o_t = int((calls_of(lk, th) == "MUT").sum())
                e_t = float(L.p_mut(lk["n"].values, th, 1, rho).sum())
                base.update({f"P_s_{arm}": st["P_s"], f"c_s_{arm}": st["c_s"], f"x_s_{arm}": st["x_s"],
                             f"n_s_{arm}": st["n_s"], f"q_T_{arm}": st["q_T"], f"flag_own_{arm}": st["flag"],
                             f"leuk_covered_{arm}": len(lk), f"leuk_umis_{arm}": int(lk["n"].sum()),
                             f"leuk_mut_umis_{arm}": int(lk["MUT"].sum()), f"neg_covered_{arm}": len(ngc),
                             f"O_T_{arm}": o_t, f"E_T_{arm}": e_t, f"FDR_{arm}": e_t / o_t if o_t else float("inf")})
                if T in C.POOL[S]:
                    for si, (A, B) in enumerate(splits):
                        a, b = ngc[ngc["bc"].isin(A)], ngc[ngc["bc"].isin(B)]
                        xA, nA = int(a["MUT"].sum()), int(a["n"].sum())
                        thA = theta_v3(st["P_s"], xA, nA)
                        pr = L.p_mut(b["n"].values, thA, 1, thA)
                        ts.append(dict(sample=S, target=T, seed=si, arm=arm, indel=T in indel_T,
                                       flagged=bool(fl["flag"]), nA=nA, B_cov=len(b),
                                       O=int((calls_of(b, thA) == "MUT").sum()), probs=pr, E=float(pr.sum())))
            trows.append(base)

    tdf = pd.DataFrame(trows)

    # ---------------- §5.1 pooled calibration, full pool, both runs under the v3 rule ----------------
    def as89(arm):  # adapt rows to score89.pooled's O_<c>/probs_<c> keys
        return [dict(r, **{"O_x": r["O"], "probs_x": r["probs"]}) for r in ts if r["arm"] == arm]

    prim = lambda r: (not r["flagged"]) and r["nA"] > 0  # noqa: E731
    frag = lambda r: (not r["flagged"]) and r["nA"] >= C.FRAG_MIN_NA  # noqa: E731
    sub_indel = lambda r: prim(r) and r["indel"]  # noqa: E731
    psum, by_seed = [], []
    for arm in ARMS:
        rows = as89(arm)
        for lab, inc in [("primary", prim), ("fragility", frag), ("indel_subpool", sub_indel)]:
            seeds = [pooled(rows, "x", si, inc) for si in range(len(C.SEEDS))]
            for si, x in enumerate(seeds):
                by_seed.append(dict(test=lab, arm=arm, seed=si, **x))
            res = {k: med([x[k] for x in seeds]) for k in ["O", "E", "p_up", "p_low", "ratio", "units"]}
            psum.append(dict(test=lab, arm=arm,
                             **{f"{k}_{q}": v for k, vv in res.items() for q, v in zip(["med", "lo", "hi"], vv)},
                             pass_=bool(res["p_up"][0] >= 0.05 and res["ratio"][0] <= 1.5),
                             seeds_pass=float(np.mean([x["p_up"] >= 0.05 and x["ratio"] <= 1.5 for x in seeds]))))
    psum = pd.DataFrame(psum)
    # sanity: the v3 arm reproduces #89's v3 pooled numbers exactly
    for lab in ("primary", "fragility"):
        a = psum[(psum["test"] == lab) & (psum["arm"] == "v3")].iloc[0]
        b = pooled89[(pooled89["test"] == lab) & (pooled89["contender"] == "v3")].iloc[0]
        for k in ("O_med", "E_med", "p_up_med", "ratio_med"):
            assert abs(a[k] - b[k]) <= 1e-9 * max(1, abs(b[k])), f"v3 arm != #89 {lab} {k}: {a[k]} vs {b[k]}"
    checks["v3_arm_reproduces_89_pooled"] = True

    # ---------------- §5.3 site screen per arm ----------------
    flat = pd.DataFrame([{k: v for k, v in r.items() if k != "probs"} for r in ts])
    site_rows = []
    for (S, T), g in flat.groupby(["sample", "target"], sort=False):
        row = dict(sample=S, target=T, indel=bool(g["indel"].iloc[0]), flagged=bool(g["flagged"].iloc[0]))
        for arm in ARMS:
            inc = g[(g["arm"] == arm) & (g["nA"] > 0)]
            row[f"included_seeds_{arm}"] = len(inc)
            ps, rr = [], []
            for _, r in inc.iterrows():
                o, e = int(r["O"]), float(r["E"])
                ps.append(0.0 if (e == 0 and o > 0) else float(poisson.sf(o - 1, e)) if o > 0 else 1.0)
                rr.append(L.ratio(o, e))
            for k, v in [("O_B", inc["O"].tolist()), ("E_B", inc["E"].tolist()), ("OE", rr), ("site_p", ps)]:
                row[f"{k}_{arm}_med"] = med(v)[0] if len(v) else np.nan
        site_rows.append(row)
    site = pd.DataFrame(site_rows)
    exc = {}
    for arm in ARMS:
        tested = site[(~site["flagged"]) & (site[f"included_seeds_{arm}"] > 0)]
        h = holm({(r["sample"], r["target"]): r[f"site_p_{arm}_med"] for _, r in tested.iterrows()}, C.ALPHA_HOLM)
        site[f"holm_adj_{arm}"] = [h.get((r["sample"], r["target"]), (np.nan, None))[0] for _, r in site.iterrows()]
        site[f"site_exception_{arm}"] = [bool(h.get((r["sample"], r["target"]), (np.nan, False))[1])
                                         for _, r in site.iterrows()]
        exc[arm] = sorted(f"{s} {t}" for (s, t), (_, rej) in h.items() if rej)
    new_indel_exc = sorted(set(site.loc[site["indel"] & site["site_exception_proto"], "sample"] + " " +
                               site.loc[site["indel"] & site["site_exception_proto"], "target"]) - set(exc["v3"]))
    assert exc["v3"] == ["AML5 T11"], f"v3 site exceptions drifted from #89: {exc['v3']}"

    # ---------------- §8 criteria ----------------
    tdf = tdf.set_index(["sample", "target"])
    gain_rows = []
    for S, T in GAIN_SET:
        r = tdf.loc[(S, T)]
        status = "flagged" if r["flagged"] else ("no-signal" if r.get("no_signal") is True else "gated")
        gain_rows.append(dict(sample=S, target=T, status=status, O_T_v3=int(r["O_T_v3"]), O_T_proto=int(r["O_T_proto"]),
                              FDR_v3=r["FDR_v3"], FDR_proto=r["FDR_proto"], band=bands[S]["band"]))
    gain = pd.DataFrame(gain_rows)
    g = gain[gain["status"] == "gated"]
    o_v, o_p = int(g["O_T_v3"].sum()), int(g["O_T_proto"].sum())
    band_sum = int(sum(bands[S]["band"] for S in g["sample"].unique()))
    rel = (o_p / o_v - 1) if o_v else float("nan")
    c1 = bool((o_p - o_v) > band_sum and rel >= 0.10)
    pp = psum[(psum["test"] == "primary") & (psum["arm"] == "proto")].iloc[0]
    c2 = bool(pp["pass_"] and not new_indel_exc)
    c3 = bool((g["FDR_proto"] <= 0.2).all()) if len(g) else False
    # sensitivity (non-gating): include sAML1 T11 despite its flag, and its band
    gs = gain
    sens = dict(O_T_v3=int(gs["O_T_v3"].sum()), O_T_proto=int(gs["O_T_proto"].sum()),
                band_sum=int(sum(bands[S]["band"] for S in gs["sample"].unique())))
    reasons = []
    if not c1:
        reasons.append(f"gain {o_p - o_v:+d} calls ({rel:+.1%}) vs summed tie bands {band_sum} and the 10 % bar")
    if not c2:
        reasons.append("calibration: " + ("prototype fails §5.1" if not pp["pass_"] else "") +
                       (f"; new indel site exception(s) {new_indel_exc}" if new_indel_exc else ""))
    if not c3:
        reasons.append("precision: FDR > 0.2 at " + ", ".join(f"{s} {t}" for s, t, f in
                                                             zip(g["sample"], g["target"], g["FDR_proto"]) if f > 0.2))
    verdict = dict(
        rule="v3 (primary verdict #90: rule 4, v3 stays provisional)",
        gain_set=gain_rows, gated=[f"{s} {t}" for s, t in zip(g["sample"], g["target"])],
        O_T_gated_v3=o_v, O_T_gated_proto=o_p, gain_calls=o_p - o_v, gain_rel=rel, band_sum=band_sum,
        criterion1_gain=c1, criterion2_calibration=c2, criterion3_precision=c3,
        proto_pooled_pass=bool(pp["pass_"]), site_exceptions=exc, new_indel_site_exceptions=new_indel_exc,
        worth_code_ticket=bool(c1 and c2 and c3), reasons=reasons, sensitivity_incl_flagged_T11=sens,
        bands=bands, checks=checks)

    tdf.reset_index().to_csv(f"{out}/targets.tsv", sep="\t", index=False)
    psum.to_csv(f"{out}/pooled.tsv", sep="\t", index=False)
    pd.DataFrame(by_seed).to_csv(f"{out}/pooled_by_seed.tsv", sep="\t", index=False)
    site.to_csv(f"{out}/site_screen.tsv", sep="\t", index=False)
    gain.to_csv(f"{out}/gain.tsv", sep="\t", index=False)
    json.dump(verdict, open(f"{out}/verdict.json", "w"), indent=1, default=str)
    pd.set_option("display.width", 250)
    print(psum.to_string())
    print(gain.to_string())
    print(json.dumps({k: v for k, v in verdict.items() if k not in ("checks",)}, indent=1, default=str))


if __name__ == "__main__":
    main()
