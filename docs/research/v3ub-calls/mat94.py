#!/usr/bin/env python3
"""Ticket #94 (map #92): materialise per-cell v3-ub calls for the 8 confirmation-pool patients.

Rule: P_null_ub = min(1, P_s + common89.cpu(n_neg_mut, n_neg_umis)) from the v3 run's <S>_targets.csv (all
negatives, no A/B split), then v3's imported `call` (genotype-v3@61d8809) with min_mut_umis = 1.
Writes OUTDIR/<S>_v3ub.csv = v3 per-cell columns + P_null_ub, p_ub, genotype_ub (every row of the v3 table).
Checks (aggregate output only; targets as T1..Tn):
  1. re-calling v3 with its own P_null reproduces `genotype` (and `p`) with 0 differences;
  2. gated tumour-call sum (scoring's QC + Leukemic sets, phased locus once as the median) = v3 2,534 / v3-ub 2,249;
  3. outside the gated set (unflagged, leuk-covered) = 775 / 601;
  plus per-target O_T agreement with q89_score/out/targets.tsv.
Usage: mat94.py OUTDIR Q89_OUT
"""
import json
import sys

import numpy as np
import pandas as pd

import common89 as C
import score_lib as L


def ndiff(a, b):
    a, b = pd.Series(list(a), dtype=object), pd.Series(list(b), dtype=object)
    return int((~((a == b) | (a.isna() & b.isna()))).sum())


def pdiff(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return int((~(np.isclose(a, b, rtol=1e-9, atol=0) | (np.isnan(a) & np.isnan(b)))).sum())


def locus_sum(df, col):
    tot = 0.0
    for _, g in df.groupby("locus", sort=False):
        tot += float(np.median(g[col])) if len(g) > 1 else float(g[col].iloc[0])
    return tot


def main():
    out, q89 = sys.argv[1], sys.argv[2]
    C.check_fixed_sets()
    tum = pd.read_csv(f"{q89}/tumour.tsv", sep="\t")
    frozen = pd.read_csv(f"{q89}/targets.tsv", sep="\t").set_index(["sample", "target"])
    gated = set(zip(tum.loc[tum["status"] == "gated", "sample"], tum.loc[tum["status"] == "gated", "target"]))

    rows, report = [], {}
    for S in C.PATIENTS:
        base = f"{C.RUNS}/v3_{S}/results/{S}/{S}"
        full = pd.read_csv(f"{base}.csv")
        tg = pd.read_csv(f"{base}_targets.csv").set_index("target_id")
        pn_ub = {t: min(1.0, float(r["P_s"]) + C.cpu(int(r["n_neg_mut"]), int(r["n_neg_umis"])))
                 for t, r in tg.iterrows()}
        full["P_null_ub"] = full["target_id"].map(pn_ub)
        mut, wt = full["MUT"].values.astype(int), full["WT"].values.astype(int)
        p_re, g_re = L.call(mut, wt, full["P_null"].values.astype(float), 1)
        full["p_ub"], full["genotype_ub"] = L.call(mut, wt, full["P_null_ub"].values.astype(float), 1)
        chk1 = dict(rows=len(full), genotype_diff=ndiff(g_re, full["genotype"]), p_diff=pdiff(p_re, full["p"]))
        full.to_csv(f"{out}/{S}_v3ub.csv", index=False)

        # scoring's view: QC cells, T labels, sets (common89.load), with the new call joined on (cell, target_id)
        d, tgT, _, _ = C.load(S)
        d = d.merge(full[["cell", "target_id", "genotype_ub"]], on=["cell", "target_id"], how="left", validate="1:1")
        neg_diff = 0
        for T in tgT.index:
            st = C.target_stats(d, tgT, T)
            r = tgT.loc[T]
            neg_diff += int((st["x_s"], st["n_s"]) != (int(r["n_neg_mut"]), int(r["n_neg_umis"])))
            lk = d[(d["T"] == T) & (d["set"] == "leuk") & (d["n"] >= 1)]
            o3, oub = int((lk["genotype"] == "MUT").sum()), int((lk["genotype_ub"] == "MUT").sum())
            fz = frozen.loc[(S, T)]
            rows.append(dict(sample=S, target=T, flagged=st["flag"], leuk_covered=len(lk), gated=(S, T) in gated,
                             O_T_v3=o3, O_T_ub=oub, frozen_v3=int(fz["O_T_v3"]), frozen_ub=int(fz["O_T_v3-ub"])))
        report[S] = dict(check1=chk1, targets=len(tg), targets_neg_counts_differ_from_scoring=neg_diff,
                         MUT_v3=int((full["genotype"] == "MUT").sum()), MUT_ub=int((full["genotype_ub"] == "MUT").sum()),
                         WT_v3=int((full["genotype"] == "WT").sum()), WT_ub=int((full["genotype_ub"] == "WT").sum()))

    t = pd.DataFrame(rows)
    t["locus"] = [f"{s} {{T3-T5}}" if (s == C.PHASED[0] and x in C.PHASED[1]) else f"{s} {x}"
                  for s, x in zip(t["sample"], t["target"])]
    g = t[t["gated"]]
    o = t[(~t["gated"]) & (~t["flagged"]) & (t["leuk_covered"] > 0)]
    res = dict(per_patient=report,
               check1_pass=all(r["check1"]["genotype_diff"] == 0 for r in report.values()),
               gated=dict(v3=locus_sum(g, "O_T_v3"), ub=locus_sum(g, "O_T_ub"), targets=len(g)),
               outside=dict(v3=int(o["O_T_v3"].sum()), ub=int(o["O_T_ub"].sum()), targets=len(o)),
               per_target_O_T_mismatch=int(((t["O_T_v3"] != t["frozen_v3"]) | (t["O_T_ub"] != t["frozen_ub"])).sum()))
    res["check2_pass"] = (res["gated"]["v3"], res["gated"]["ub"]) == (2534, 2249)
    res["check3_pass"] = (res["outside"]["v3"], res["outside"]["ub"]) == (775, 601)
    t.drop(columns=["locus"]).to_csv(f"{out}/check_targets.tsv", sep="\t", index=False)
    json.dump(res, open(f"{out}/checks.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
