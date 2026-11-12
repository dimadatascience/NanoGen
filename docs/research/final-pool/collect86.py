#!/usr/bin/env python3
"""Ticket 86: final pool for the pooled calibration test (aggregate only, targets T1..Tn).
Imports #76's collect.py and power_check.py unchanged (q76_power/scripts, = 9a2d94f).
usage: collect86.py <q86 out dir>   (expects raw_/wc_/cen_<S>.tsv for all eight patients)"""
import sys

import pandas as pd

import collect
import power_check as pc

d = sys.argv[1]
CENSUS = ["PALM07_E", "PALM08_E", "PALM11_E", "PALM14_E", "PALM21_E"]
RESTART = ["AML4", "AML5", "sAML1"]
ALL = CENSUS + RESTART
GRID = pc.C_GRID
C1 = "E_B_c0.01_cen"
STOP_R80 = 2.0
SCREEN_E = 5.7

print("## Per-target table (collect.py census, unchanged; central scale for E_B)\n")
collect.census(d, ALL)
x = pd.read_csv(f"{d}/census.tsv", sep="\t")
cols = [f"E_B_c{c}_cen" for c in GRID]
# E_B is linear-ish in coverage: >0 at one c iff >0 at all c; check it.
pos = x[cols] > 0
assert (pos.all(axis=1) == pos.any(axis=1)).all(), "E_B > 0 not consistent across the c grid"
x["in_pool"] = pos.all(axis=1)

print("\n## Verdicts (descriptive only; the pooled test replaces per-target evaluability)\n")
print(x.groupby(["sample", "verdict"]).size().unstack(fill_value=0).reindex(ALL).fillna(0).astype(int).to_string())

print("\n## Pool membership: targets with E_B > 0\n")
mem = x[x["in_pool"]].groupby("sample", sort=False)["target"].apply(lambda t: ", ".join(t)).reindex(ALL)
out_ = x[~x["in_pool"]].groupby("sample", sort=False)["target"].apply(lambda t: ", ".join(t)).reindex(ALL)
print(collect.md(pd.DataFrame({"sample": ALL,
                               "targets_total": [int((x["sample"] == s).sum()) for s in ALL],
                               "in_pool": [int(((x["sample"] == s) & x["in_pool"]).sum()) for s in ALL],
                               "in_pool_targets": mem.fillna("-").values,
                               "E_B_0_targets": out_.fillna("-").values})))

print("\n## Per-target E_B (central scale) across the c grid, pool members\n")
y = x[x["in_pool"]][["sample", "target", "neg_covered", "B_cells_med_cen"] + cols + ["verdict"]]
print(collect.md(y))
scr = x[x[C1] >= SCREEN_E]
print(f"\nTargets with E_B >= {SCREEN_E} at c = 1 % (screen's 80 % power point for a 3x site miss): {len(scr)}"
      + ("" if scr.empty else " -> " + ", ".join(f"{a} {b} ({e:.2f})" for a, b, e in scr[["sample", "target", C1]].itertuples(index=False))))


def pool(df, label):
    rows = []
    for c in GRID:
        E = float(df[f"E_B_c{c}_cen"].sum())
        p = pc.power_numbers(E)
        rows.append(dict(pool=label, c=c, E_B=round(E, 2), r80=p["r80"]))
    return rows


rows = []
for s in ALL:
    rows += pool(x[x["sample"] == s], s)
rows += pool(x[x["sample"].isin(CENSUS)], "census (5)")
rows += pool(x[x["sample"].isin(RESTART)], "restart (3)")
rows += pool(x, "FINAL POOL (8)")
for s in ALL:
    rows += pool(x[x["sample"] != s], f"final pool minus {s}")
t = pd.DataFrame(rows)
order = list(dict.fromkeys(t["pool"]))


def wide(v):
    w = t.pivot(index="pool", columns="c", values=v).reindex(order).reset_index()
    return w.rename(columns=lambda k: f"{v} c={k}" if k != "pool" else k)


print("\n## Pools (E_B summed over pool targets, central scale; r80 from power_check.power_numbers)\n")
print(collect.md(wide("E_B")))
print()
print(collect.md(wide("r80")))
share = (x[x["in_pool"]].groupby("sample")[C1].sum() / x[C1].sum()).reindex(ALL)
print("\nShare of pooled E_B at c = 1 %: " + ", ".join(f"{s} {v:.1%}" for s, v in share.items()))
fin = t[(t["pool"] == "FINAL POOL (8)") & (t["c"] == 0.01)].iloc[0]
print(f"\nE* = {pc.e_star():.1f}")
print(f"\n## Stop rule\n\nFinal pool r80 at c = 1 % = {fin['r80']:.2f} (E_B {fin['E_B']:.2f}); stop if > {STOP_R80}: "
      + ("**STOP**" if fin["r80"] > STOP_R80 else "**PASS (do not stop)**"))
lo = t[(t["pool"].str.startswith("final pool minus")) & (t["c"] == 0.01)]
print(f"Leave-one-patient-out r80 at c = 1 %: max {lo['r80'].max():.2f} ({lo.loc[lo['r80'].idxmax(), 'pool']}); "
      f"all <= {STOP_R80}: {bool((lo['r80'] <= STOP_R80).all())}")
x.to_csv(f"{d}/pool86.tsv", sep="\t", index=False)
t.to_csv(f"{d}/pools86.tsv", sep="\t", index=False)
