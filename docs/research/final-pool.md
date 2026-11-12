# Final pool for the pooled calibration test (ticket 86, map 72)

**Answer.** The final pool is **8 patients and 57 targets**. Pooled E_B at c = 1 % is **38.0**, and **r80 = 1.45**. The stop rule (stop if r80 > 2.0 at c = 1 %) **does not trigger**. It still holds with any single patient left out: the worst case is r80 1.55, without PALM21E. No target reaches the screen's 3× power point (E_B ≥ 5.7 at c = 1 %).

Aggregate numbers only, targets as T1…Tn (key in the HPC README). HPC folder `q86_pool/`, job 9074545 (an earlier run, 9074539, stopped at a byte-identity check; see Method).

## Method (fixed by the ticket)

- **Proxy:** `power_check.py` at `9a2d94f`, unchanged (md5 `1dc67c18…` on the HPC, equal to `git show 9a2d94f:docs/research/power-check/power_check.py`). Scenario grid c = 0.25 / 0.5 / 1 / 2 / 4 %. E_B is at the central proxy→v3 scale (0.9073 UMIs, 0.9513 cells), as in #76 and #84. Verdicts use raw and worst-case (0.685 / 0.709) scales and the unchanged `collect.py` rule.
- **Census inputs: final canonical BAMs.** PALM08E, 11E, 14E and 21E use the BAMs regenerated in ticket 85. PALM07E uses library 2 (ticket 88), with its derived library-2 QC list, negatives and labels. BEDs, QC lists, negatives and labels are as in the HPC README.
- **Re-run:** per-cell coverage on these BAMs had already been computed by the same code in q85 and q88. q86 re-summarised it at all three scales.
  - The re-summarised raw table must equal the q85 / q88 raw table: string fields identical, numbers within 1e-9 (`same86.py`).
  - All five match, with a maximum absolute difference of 8.9e-16 (last float bit, from summation order). The first submission, which demanded byte identity, stopped on exactly that bit.
- **Restart patients:** AML4, AML5 and sAML1 passed intake (ticket 84). Their tables are taken unchanged from `q84_restart/out/full`.
- **Pool:** sum of E_B over every target with E_B > 0. E_B > 0 holds at every c or at none, which `collect86.py` asserts. r80 comes from `power_check.power_numbers`.

## Pool membership

| patient | targets | in pool (E_B > 0) | E_B = 0 (outside the pool) |
|---|---|---|---|
| PALM07E | 5 | T1–T5 | – |
| PALM08E | 11 | T1, T3–T7 | T2, T8–T11 |
| PALM11E | 5 | T1, T3–T5 | T2 |
| PALM14E | 7 | T1–T7 | – |
| PALM21E | 12 | T1–T5, T7–T12 | T6 |
| AML4 | 11 | T1–T4, T6, T9, T10 | T5, T7, T8, T11 |
| AML5 | 16 | T1–T4, T7, T10, T11, T16 | T5, T6, T8, T9, T12–T15 |
| sAML1 | 16 | T1, T2, T4, T7–T12 | T3, T5, T6, T13–T16 |
| **total** | **83** | **57** | 26 |

PALM05E is not in the pool (no annotation; see the map's fog).

## Pooled E_B and r80

| pool | E_B 0.25 % | 0.5 % | **1 %** | 2 % | 4 % | r80 0.25 % | 0.5 % | **1 %** | 2 % | 4 % |
|---|---|---|---|---|---|---|---|---|---|---|
| PALM07E | 1.67 | 3.28 | 6.28 | 10.7 | 15.3 | 4.02 | 2.78 | 2.36 | 1.91 | 1.77 |
| PALM08E | 1.31 | 2.49 | 4.37 | 6.43 | 7.68 | 4.22 | 3.19 | 2.61 | 2.30 | 2.07 |
| PALM11E | 1.20 | 2.25 | 4.09 | 6.36 | 6.51 | 4.61 | 3.52 | 2.79 | 2.33 | 2.28 |
| PALM14E | 0.69 | 1.34 | 2.40 | 4.05 | 4.71 | 6.24 | 4.12 | 3.30 | 2.82 | 2.66 |
| PALM21E | 3.25 | 6.23 | 11.6 | 20.8 | 35.1 | 2.80 | 2.38 | 1.86 | 1.67 | 1.47 |
| AML4 | 0.66 | 1.28 | 2.51 | 4.51 | 6.56 | 6.53 | 4.30 | 3.15 | 2.53 | 2.26 |
| AML5 | 0.59 | 1.16 | 2.28 | 3.94 | 5.69 | 7.25 | 4.78 | 3.47 | 2.60 | 2.40 |
| sAML1 | 1.23 | 2.37 | 4.49 | 8.20 | 12.0 | 4.49 | 3.34 | 2.54 | 2.08 | 1.88 |
| census (5) | 8.12 | 15.6 | 28.7 | 48.3 | 69.3 | 2.10 | 1.73 | 1.54 | 1.40 | 1.33 |
| restart (3) | 2.48 | 4.81 | 9.28 | 16.6 | 24.2 | 3.20 | 2.61 | 2.08 | 1.75 | 1.60 |
| **final pool (8)** | **10.6** | **20.4** | **38.0** | **64.9** | **93.6** | **1.93** | **1.64** | **1.45** | **1.34** | **1.28** |

- E\* = 30.9 (1.5× miss at 80 % power). The final pool passes E\* at c = 1 % but not at c ≤ 0.5 %.
- **Share of pooled E_B at c = 1 %:** PALM21E 30 %, PALM07E 17 %, sAML1 12 %, PALM08E 12 %, PALM11E 11 %, AML4 7 %, PALM14E 6 %, AML5 6 %. No single patient dominates.
- **Leave one patient out (r80 at c = 1 %):** 1.48–1.55. The worst case, without PALM21E, is E_B 26.4 and r80 1.55.
- **Changes against the #84 preview** (census on its earlier canonical BAMs: E_B 31.7, r80 1.50): the census now gives 28.7 (was 22.4). The gains come from the known-list BAMs (ticket 85) and PALM07E library 2 (ticket 88). PALM07E is 6.28 at the central scale; ticket 88's 6.92 is the raw-scale sum of the same table.

## Per-target status (descriptive only)

The pooled test replaces per-target evaluability. This table is reported, not gated.

| patient | evaluable | borderline | not evaluable |
|---|---|---|---|
| PALM07E | 1 (T5) | 1 (T1) | 3 |
| PALM08E | 0 | 0 | 11 |
| PALM11E | 0 | 1 (T1) | 4 |
| PALM14E | 0 | 0 | 7 |
| PALM21E | 0 | 2 (T8, T11) | 10 |
| AML4 | 0 | 0 | 11 |
| AML5 | 0 | 0 | 16 |
| sAML1 | 0 | 1 (T9) | 15 |

- **Largest per-target E_B at c = 1 %:** PALM07E T5 3.42; PALM11E T1 3.06; PALM21E T11 2.60; PALM07E T1 2.00.
- **No target has E_B ≥ 5.7 at c = 1 %**, so the site screen has < 80 % power for a 3× miss at every site. At c = 4 %, PALM07E T5 (6.73), T1 (5.94), PALM21E T8 (6.32) and T11 (9.70) pass 5.7. The full per-target grid is in `final-pool/out/report.md`.

## Stop rule

**Final pool r80 at c = 1 % = 1.45 ≤ 2.0: do not stop.** The map's caveats still hold:
- the pooled test can see a ≈ 1.45× average miss at 80 % power when c = 1 %, and a 2× average miss comfortably;
- neither the pool nor the site screen can rule out a T3-sized (≈ 1.7×) single-site miss.

## Files

- `final-pool/collect86.py`: pooling and report (imports `collect.py` and `power_check.py` at `9a2d94f`, unchanged).
- `final-pool/same86.py`: raw-table equality check against q85 / q88.
- `final-pool/out/report.md`, `pool86.tsv` (per target, all scales), `pools86.tsv` (pools × c).
