# Exploratory caller-only sweep (PALM22E, PALM09E)

> **EXPLORATORY — outside the pre-registered verdict.** Wayfinder ticket *Exploratory caller-only sweep* (#68), map *T/B negative controls vs blast calls* (#60). The computations are those of the pre-registration (`docs/research/negative-control-preregistration.md` @ `76c3f99`). Its decision rules are applied here for orientation only. The verdict ticket (#69) reads the scoring ticket (#67), not this file.

**Patient data:** aggregate numbers only, targets as T1…Tn.

## Question

Rescore the genotype-v3 run's consensus counts (#66) with caller variants outside the verdict, and report the same metrics as the scoring ticket. The variants are:

- the upper 95 % bound of `c_s`;
- `c_s` from T/B only vs T/B/NK;
- v3 + `min_mut_umis = 2`.

## Answer in brief

| variant | what changes vs v3 | effect |
|---|---|---|
| **v3-ub** (θ = P_s + one-sided 95 % Clopper–Pearson upper bound of c) | larger P_null | **Best-calibrated.** Calibration passes in 20/20 seeds at PALM22E T4, T5 and PALM09E T1, and in 17/20 at PALM09E T3; v3 passes in 5–17/20. It fixes PALM09E T3, the one unflagged target where v3 fails calibration. It costs ≤ 3.5 % of tumour MUT calls (PALM22E T5 0, T3 −6, T4 −42, T7 −13). |
| **v3-tb** (c from T/B only, NK dropped) | c almost unchanged | **No gain.** NK contamination is in line with T cells. Where negatives are heterogeneous, the high-c subgroup is B cells, not NK (PALM22E T3 B 4.9 %, PALM09E T3 B 13.5 % vs T 1.4 % / 3.4 %). Calls are identical to v3 within ±1. |
| **v3-min2** (m_min = 2) | single-UMI MUT calls removed | **Collapses to v2 + min2 in the tumour compartment.** With P_null ≈ 0.5 %, m ≥ 2 dominates the rule. Tumour MUT calls drop by more than half (PALM22E T5 2,671 → 1,214; T4 2,514 → 1,217; detection efficiency 0.74 → 0.34). Criterion 3a is a tie by construction. |

**Cross-cutting finding.** PALM09E T1 is the loss flagged in #66 (v3: 34 MUT calls vs `main` 1,272). It is **not** a call-rule effect. In the v3 run's counts, the pooled MUT fraction over covered Leukemic cells is q_T = 0.42 %, close to c_s = 0.09 %. Every variant with m_min = 1 calls the same 34 cells, and no rule over these counts can recover the clone. The loss sits upstream, in consensus/indel labelling. This answers the map's fog item *Why v3 loses MUT calls at PALM09E's indel targets* as caller-side, not statistical.

## Contenders

All contenders use the v3 run's counts, QC cells only, and P_s from the run.

| name | nominal null θ | m_min | role |
|---|---|---|---|
| v2 | P_s | 1 | reference |
| v2 + min2 | P_s | 2 | reference; 3a baseline |
| v3 | P_s + c | 1 | reference (= the v3 run) |
| v3-ub | P_s + CP₉₅(c) | 1 | sweep |
| v3-tb | P_s + c over T ∪ B (`aggregated_ct` ≠ `NK_cells`) | 1 | sweep |
| v3-min2 | P_s + c | 2 | sweep |

c is c_A in the held-out split and c_s in-sample. `main` is not part of this sweep: its counts differ, and it is scored in #67.

**Choices made here, which the pre-registration does not fix:**

- **"Upper 95 % bound"** is the one-sided 95 % Clopper–Pearson bound, `beta.ppf(0.95, x+1, n−x)`. It is exact and defined at x = 0. With n = 0 it falls back to 0, as v3 does.
- **v3-tb** uses the same A/B split of all T/B/NK negatives. It estimates c_A from A minus its NK cells, and is scored on all of B, NK included, so its calibration cells match the other contenders'.
- **Clonal targets** for criteria 2 and 3 are **provisional**. *Clonal targets in PALM22E and PALM09E* (#63) is still open, so this sweep uses its default: PALM22E T5, T3, T7; PALM09E T1, T4.
- The **one-UMI test** (§6) does not depend on the call rule, so it is left to #67.
- **Estimated true null** (FDR, efficiency) is P_s + c_s from all negatives for every contender, as pre-registered. v3-ub and v3-tb change only the rule, not the estimate of truth.

**Implementation checks** (all passed; `sweep-out/checks_*.json`):

- `call` is imported from `genotype-v3`@`61d8809` (md5-identical file on the HPC).
- Re-calling every row with the run's P_null reproduces the `genotype` column exactly.
- c_s recomputed over QC negatives matches the run to < 1e-16, and re-calling QC rows with it reproduces the calls exactly.
- All negatives are QC cells. NK cells (366 in PALM22E, 71 in PALM09E) all lie inside the T/B/NK list, and every negative is T, B or NK.
- A synthetic dataset ran end to end first.
- Software: numpy 2.0.1, in the pinned v3 image (`ghcr.io/dimadatascience/nanogen@sha256:764368d6…`). Slurm job 9002897 on `short`, 1 CPU / 4 GB.

## Results

The full per-target × contender tables, every §6 column included, are `negative-control/sweep-out/sweep_<S>.tsv`. The contamination-by-subgroup tables are `c_by_group_<S>.tsv`. Held-out values are seed medians over evaluable seeds.

### Evaluability and flags

| sample | target | evaluable seeds | c_s | q_T | flag | clonal (provisional) |
|---|---|---|---|---|---|---|
| PALM22E | T3 | 6 / 20 | 0.018 | 0.52 | – | yes |
| PALM22E | T4 | 20 | 0.006 | 0.17 | – | no (excluded pending #63) |
| PALM22E | T5 | 20 | 0.005 | 0.49 | – | yes |
| PALM22E | T7 | 0 | 0.026 | 0.51 | – | yes |
| PALM09E | T1 | 20 | 0.0009 | 0.004 | – | yes |
| PALM09E | T2 | 20 | 0.0034 | 0.0034 | **c_s ≥ 0.5 q_T** | no |
| PALM09E | T3 | 20 | 0.040 | 0.13 | – (0.040 < 0.065) | no |
| PALM09E | T4 | 0 | 0.007 | 0.09 | – | yes |

PALM22E T1, T2 and T6 have ≤ 10 covered negatives and are not evaluable. PALM22E T3 (borderline in #66) falls short in 14 of 20 seeds, so it is **not evaluable**. PALM09E T3, the germline/artefact suspect, is **not** flagged by the pre-registered rule. PALM09E T2 is.

### Criterion 1: held-out calibration (seed median)

| sample | target | v2 | v2 + min2 | v3 | **v3-ub** | v3-tb | v3-min2 |
|---|---|---|---|---|---|---|---|
| PALM22E | T4 | O/E 28.8, p 0 ✗ | 0 obs ✓ (0.55) | 0.91, 0.71 ✓ (0.85) | **0.69, 0.98 ✓ (1.0)** | 0.94, 0.65 ✓ (0.85) | 0 obs ✓ (0.55) |
| PALM22E | T5 | 47.9, 1e-6 ✗ | 0 obs ✓ (1.0) | 0.97, 0.59 ✓ (0.7) | **0.45, 0.98 ✓ (1.0)** | 1.25, 0.37 ✓ (0.6) | 0 obs ✓ (1.0) |
| PALM09E | T1 | 5.85, 0.12 ✗ | 0 obs ✓ (1.0) | 0.98, 0.64 ✓ (0.55) | **0.27, 0.98 ✓ (1.0)** | 0.96, 0.65 ✓ (0.55) | 0 obs ✓ (1.0) |
| PALM09E | T3 | 402, 3e-16 ✗ | 2.4e4, 4e-5 ✗ (0.35) | **1.70, 0.048 ✗ (0.25)** | **1.29, 0.27 ✓ (0.85)** | 1.69, 0.063 ✗ (0.25) | 1.46, 0.50 ✓ (0.55) |
| PALM09E | T2 (flagged) | 34.8 ✗ | ✓ | 1.03 ✓ | 0.37 ✓ | 0.95 ✓ | ✓ |

Each cell reads median O_B/E_B, median p_cal, then pass ✓ / fail ✗, with the fraction of seeds passing in brackets. "0 obs" means O_B = 0 MUT calls, so ratio 1 and p 1.

- At PALM09E T3, v3's pooled c (4.0 %) understates the B-cell contamination (13.5 %, 5/37 UMIs), so the tail of false MUT calls is over-dispersed. The upper bound absorbs it; dropping NK does not.
- v3-ub over-corrects slightly: O/E is 0.27–0.69 at the deep targets.

### Criterion 2 and sensitivity: tumour compartment, in-sample

| sample | target | metric | v2 | v2 + min2 | v3 | v3-ub | v3-tb | v3-min2 |
|---|---|---|---|---|---|---|---|---|
| PALM22E | T5 (clonal) | O_T | 2,671 | 1,214 | 2,671 | 2,671 | 2,671 | 1,214 |
| | | FDR | 0.017 | 0.0002 | 0.017 | 0.017 | 0.017 | 0.0002 |
| | | efficiency (π̂ 1.00, f̂ 0.49) | 0.74 | 0.34 | 0.74 | 0.74 | 0.74 | 0.34 |
| PALM22E | T3 (clonal) | O_T / FDR / eff. | 949 / 0.044 / 0.63 | 185 / 0.002 / 0.12 | 948 / 0.043 / 0.63 | 942 / 0.040 / 0.63 | 947 / 0.042 / 0.63 | 185 / 0.002 / 0.12 |
| PALM22E | T7 (clonal) | O_T / FDR / eff. | 375 / 0.10 / 0.98 | 195 / 0.005 / 0.54 | 371 / 0.083 / 0.98 | 358 / 0.052 / 0.96 | 371 / 0.083 / 0.98 | 195 / 0.005 / 0.54 |
| PALM22E | T4 | O_T / FDR / eff. | 2,530 / 0.065 / 0.61 | 1,217 / 0.004 / 0.29 | 2,514 / 0.057 / 0.60 | 2,472 / 0.050 / 0.59 | 2,514 / 0.057 / 0.60 | 1,217 / 0.004 / 0.29 |
| PALM09E | T1 (clonal) | O_T / FDR | 34 / 0.25 | 0 / ∞ | 34 / 0.25 | 34 / 0.25 | 34 / 0.25 | 0 / ∞ |
| PALM09E | T4 (clonal) | O_T / FDR | 17 / 0.13 | 5 / 0.004 | 17 / 0.13 | 16 / 0.09 | 17 / 0.13 | 5 / 0.004 |
| PALM09E | T3 | O_T / FDR | 13 / 1.05 | 2 / 1.28 | 13 / 0.75 | 12 / 0.63 | 13 / 0.75 | 2 / 0.17 |

v3 − (v2 + min2), criterion 3a (≥ 40 calls):

- PALM22E: T5 +1,457, T3 +763, T7 +176 for v3 (all pass), and similar for v3-ub (+1,457 / +757 / +163) and v3-tb.
- PALM09E: T1 +34 and T4 +12 for v3 (both **ties**).
- v3-min2: 0 everywhere (a tie by construction).

Undetermined fraction in Leukemic cells (3c): ≈ 0.10–0.34 at m_min = 1, rising to 0.41–0.83 at m_min = 2. Wilson CIs are in the TSVs.

PALM09E T1's efficiency is not evaluable (f̂ = 0.004 < ρ₀ + 0.05), and π̂ is meaningless there. The mixture finds no MUT component in these counts.

### Exploratory sample verdicts (pre-registered rule, provisional clonal set)

| contender | PALM22E | PALM09E |
|---|---|---|
| v2 | fail (1: T4, T5) | fail (1, 2) |
| v2 + min2 | pass | fail (1: T3; 2: T1 FDR ∞) |
| v3 | **pass** | fail (1: T3; 2: T1; 3a: T1, T4 ties) |
| **v3-ub** | **pass** | fail (2: T1; 3a: T1, T4 ties); **criterion 1 passes** |
| v3-tb | pass | fail (as v3) |
| v3-min2 | fail (3a ties) | fail (2: T1; 3a ties); criterion 1 passes |

PALM09E fails every contender for the same reason: T1 carries almost no MUT signal in v3's counts, so criteria 2 and 3a cannot be met there by **any** call rule.

## What this means for the map

- **v3-ub is the strongest alternative.** It strictly improves specificity, including the one calibration failure of v3 on an unflagged target, at ≤ 3.5 % sensitivity cost. Worth carrying into the verdict ticket's discussion as a candidate setting. It is outside the pre-registered contenders, so adopting it would need its own confirmation.
- **T/B only is not worth pursuing.** NK cells are not the contamination outlier. B-cell heterogeneity at the low-depth targets is the over-dispersion source.
- **min_mut_umis = 2 on top of v3 is dominated.** It halves tumour sensitivity, and v3's c already removes the vagrant single-UMI calls in negatives.
- **PALM09E T1/T2 loss is caller-side (counts), not statistical.** That fog item can be ticketed as a counts/indel-labelling investigation, comparing `main`'s and v3's consensus at those targets. Until it is resolved, PALM09E cannot pass under any call rule, and whether PALM09E T1 stays clonal (#63) matters for the verdict.

## Reproduce

The scripts are `negative-control/score_lib.py` (pre-registered computations, reusable by #67) and `negative-control/sweep.py`. They run in the pinned v3 image with `GENOTYPING_DIR` pointing at `bin/long_reads` of `genotype-v3`@`61d8809`. Input paths and the T-label key live only on the HPC, in the working folder's `q68_sweep/` (`scripts/run.sbatch`, `keys/`).
