#!/usr/bin/env python3
"""Ticket 84: collect the old-panel power check (aggregate only, T1..Tn) and pool E_B.
Imports #76's collect.py and power_check.py unchanged (q76_power/scripts, = 9a2d94f).
usage: collect84.py <q84 out dir> <q76 out dir>"""
import sys
import io
import contextlib

import numpy as np
import pandas as pd

import collect
import power_check as pc

d84, d76 = sys.argv[1], sys.argv[2]
NEW = ["AML4", "AML5", "sAML1"]
CENSUS = ["PALM07_E", "PALM08_E", "PALM11_E", "PALM14_E", "PALM21_E"]
GRID = pc.C_GRID

print("## Per-target table (collect.py census, unchanged)\n")
collect.census(d84, NEW)
x = pd.read_csv(f"{d84}/census.tsv", sep="\t")
print("\n## Verdicts\n")
print(x.groupby(["sample", "verdict"]).size().unstack(fill_value=0).to_string())

print("\n## Per-target E_B (central scale) across the c grid, targets with E_B > 0\n")
cols = [f"E_B_c{c}_cen" for c in GRID]
y = x[x[cols].sum(axis=1) > 0][["sample", "target", "neg_covered", "B_cells_med_cen"] + cols + ["verdict"]]
print(collect.md(y))
print("\nTargets with E_B >= 5.7 at c = 1 % (screen's 80 % power point for a 3x site miss):",
      int((x["E_B_c0.01_cen"] >= 5.7).sum()))


def pool(df, label):
    rows = []
    for c in GRID:
        E = float(df[f"E_B_c{c}_cen"].sum())
        p = pc.power_numbers(E)
        rows.append(dict(pool=label, c=c, E_B=round(E, 2), r80=p["r80"]))
    return rows


c76 = pd.read_csv(f"{d76}/census.tsv", sep="\t")
out = []
for s in NEW:
    out += pool(x[x["sample"] == s], s)
out += pool(x, "restart (3)")
out += pool(c76[c76["sample"].isin(CENSUS)], "census (5, canonical BAMs; #76)")
out += pool(pd.concat([c76[c76["sample"].isin(CENSUS)], x]), "census + restart (8), preview")
t = pd.DataFrame(out)
print("\n## Pools (E_B summed over targets, central scale; r80 from power_check.power_numbers)\n")
print(collect.md(t.pivot(index="pool", columns="c", values="E_B").reset_index().rename(columns=lambda k: f"E_B c={k}" if k != "pool" else k)))
print()
print(collect.md(t.pivot(index="pool", columns="c", values="r80").reset_index().rename(columns=lambda k: f"r80 c={k}" if k != "pool" else k)))
print(f"\nE* = {pc.e_star():.1f}")
