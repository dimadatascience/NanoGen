#!/usr/bin/env python3
"""Collect the aggregate power-check tables (no barcodes) into Markdown.

usage: collect.py calib OUT_DIR             proxy-vs-v3 gap on the calibration samples; prints scale factors
       collect.py census OUT_DIR SAMPLES... per-target evaluable / underpowered table
"""
import sys

import pandas as pd

C_REF = 0.01  # reference scenario reported in the main table (the full grid is in the TSVs)


def md(df):
    """Markdown table without the tabulate dependency."""
    def f(v):
        if isinstance(v, float):
            return f"{v:.3g}"
        return str(v)
    lines = ["| " + " | ".join(df.columns) + " |", "|" + "---|" * len(df.columns)]
    lines += ["| " + " | ".join(f(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join(lines)


def calib(d):
    x = pd.concat([pd.read_csv(f"{d}/raw_{s}.tsv", sep="\t") for s in ("PALM22_E", "PALM09_E")])
    cols = ["sample", "target", "neg_covered", "v3_neg_covered", "neg_umis", "v3_neg_umis",
            "leuk_covered", "v3_leuk_covered", "cells_exact", "cells_v3_gt_proxy",
            "eval_seeds", "eval_seeds_v3", "A_umis_med", "A_umis_med_v3", "B_cells_med", "B_cells_med_v3"]
    print(md(x[cols]))
    big = x[x["v3_neg_covered"] >= 100]
    su = big["v3_neg_umis"].sum() / big["neg_umis"].sum()
    sc = big["v3_neg_covered"].sum() / big["neg_covered"].sum()
    rng_u = (big["v3_neg_umis"] / big["neg_umis"]).agg(["min", "max"]).tolist()
    rng_c = (big["v3_neg_covered"] / big["neg_covered"]).agg(["min", "max"]).tolist()
    agree = int((x["evaluable"] == x["evaluable_v3"]).sum())
    print(f"\nscale_umis={su:.4f} (per-target {rng_u[0]:.3f}-{rng_u[1]:.3f}); "
          f"scale_cells={sc:.4f} (per-target {rng_c[0]:.3f}-{rng_c[1]:.3f}); "
          f"targets with >=100 v3-covered negatives: {len(big)}; evaluable proxy==v3: {agree}/{len(x)}")


def verdict(r):
    """The proxy is an upper bound on v3's n (no cell has v3 n > proxy on PALM22E/09E), so a raw
    failure is final; a pass that survives the worst measured proxy->v3 scale is evaluable."""
    if not r["evaluable"]:
        return "not evaluable"
    if not r["evaluable_wc"]:
        return "borderline"
    return "evaluable"


def census(d, samples):
    out = []
    for s in samples:
        x = pd.read_csv(f"{d}/raw_{s}.tsv", sep="\t")
        wc = pd.read_csv(f"{d}/wc_{s}.tsv", sep="\t")
        cen = pd.read_csv(f"{d}/cen_{s}.tsv", sep="\t")
        x["evaluable_wc"] = wc["evaluable_cal"].values
        x["eval_seeds_wc"] = wc["eval_seeds_cal"].values
        for c in ("A_umis_med_cal", "B_cells_med_cal", "leuk_covered_cal", "c_min_1.5x"):
            x[c.replace("_cal", "") + "_cen"] = cen[c].values
        for c in (0.0025, 0.005, 0.01, 0.02, 0.04):
            x[f"E_B_c{c}_cen"] = cen[f"E_B_c{c}_cal"].values
            for k in ("r80", "fail_r1", "fail_r15", "r80_lo"):
                x[f"{k}_c{c}_cen"] = cen[f"{k}_c{c}"].values
        x["verdict"] = x.apply(verdict, axis=1)
        out.append(x)
    x = pd.concat(out)
    cols = ["sample", "target", "neg_qc", "leuk_qc", "neg_covered", "B_cells_med", "A_umis_med", "eval_seeds",
            "eval_seeds_wc", "leuk_covered", "B_cells_med_cen", "A_umis_med_cen", "leuk_covered_cen",
            f"E_B_c{C_REF}_cen", f"r80_c{C_REF}_cen", f"fail_r1_c{C_REF}_cen", f"fail_r15_c{C_REF}_cen",
            "c_min_1.5x_cen", "verdict"]
    print(md(x[cols]))
    x.to_csv(f"{d}/census.tsv", sep="\t", index=False)


if __name__ == "__main__":
    {"calib": lambda: calib(sys.argv[2]), "census": lambda: census(sys.argv[2], sys.argv[3:])}[sys.argv[1]]()
