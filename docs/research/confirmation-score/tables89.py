#!/usr/bin/env python3
"""Markdown tables for docs/research/confirmation-score.md from out/*.tsv (stdlib only). Usage: tables89.py OUTDIR"""
import csv
import sys

O = sys.argv[1]


def rd(f):
    return list(csv.DictReader(open(f"{O}/{f}"), delimiter="\t"))


def g(x, p=3):
    if x in ("", "nan", None):
        return "–"
    if x == "inf":
        return "∞"
    return f"{float(x):.{p}g}"


print("POOLED")
for r in rd("pooled.tsv"):
    print(f"| {r['test']} | {r['contender']} | {g(r['O_med'])} ({g(r['O_lo'])}–{g(r['O_hi'])}) "
          f"| {g(r['E_med'])} ({g(r['E_lo'])}–{g(r['E_hi'])}) | {g(r['ratio_med'])} ({g(r['ratio_lo'])}–{g(r['ratio_hi'])}) "
          f"| {g(r['p_up_med'], 2)} ({g(r['p_up_lo'], 2)}–{g(r['p_up_hi'], 2)}) | {g(r['p_low_med'], 2)} "
          f"| {g(r['units_med'])} | {r['pass_']} | {g(r['seeds_pass'])} |")
print("PATIENTS")
P = {}
for r in rd("patients.tsv"):
    P.setdefault(r["sample"], {})[r["contender"]] = r
for s, d in P.items():
    a, b = d["v3"], d["v3-ub"]
    print(f"| {s} | {g(a['O_med'])} / {g(a['E_med'])} | {g(a['OE_med'])} ({g(a['garwood90_lo_med'])}–{g(a['garwood90_hi_med'])}) "
          f"| {g(b['O_med'])} / {g(b['E_med'])} | {g(b['OE_med'])} ({g(b['garwood90_lo_med'])}–{g(b['garwood90_hi_med'])}) "
          f"| {g(a['lopo_OE_med'])}, p {g(a['lopo_p_med'], 2)} | {g(b['lopo_OE_med'])}, p {g(b['lopo_p_med'], 2)} |")
print("TARGETS")
T = {(r["sample"], r["target"]): r for r in rd("targets.tsv")}
for r in rd("site_screen.tsv"):
    t = T[(r["sample"], r["target"])]
    st = "flagged" if r["flagged"] == "True" else (
        "**site exception (v3)**" if r["site_exception_v3"] == "True" else "tested")
    print(f"| {r['sample']} | {r['target']} | {t['neg_covered']} | {g(t['B_cov_med'])} | {g(t['A_umis_med'])} "
          f"| {g(t['c_s'], 2)} | {g(t['cpu_s'], 2)} | {g(t['c_A_med'], 2)} | {g(t['q_T'], 2)} | {r['included_seeds']} "
          f"| {g(r['O_B_v3_med'])} / {g(r['E_B_v3_med'], 2)} | {g(r['site_p_v3_med'], 2)} "
          f"| {g(r['O_B_v3-ub_med'])} / {g(r['E_B_v3-ub_med'], 2)} | {g(r['site_p_v3-ub_med'], 2)} | {st} |")
print("TUMOUR")
for r in rd("tumour.tsv"):
    print(f"| {r['sample']} | {r['target']} | {r['status']} | {r['leuk_covered']} | {r['n_s']} "
          f"| {g(r.get('pi_hat', ''), 2)} / {g(r.get('f_hat', ''), 2)} | {r['O_T_v3']} | {g(r['FDR_v3'], 2)} "
          f"| {r['O_T_v3-ub']} | {g(r['FDR_v3-ub'], 2)} | {g(r['loss_target'], 2)} "
          f"| {g(r.get('efficiency_v3', ''), 2)} / {g(r.get('efficiency_v3-ub', ''), 2)} "
          f"| {g(r['undet_leuk_v3'], 2)} / {g(r['undet_leuk_v3-ub'], 2)} |")
