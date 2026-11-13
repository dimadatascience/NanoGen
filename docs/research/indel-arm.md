# Indel arm: is position-tolerant matching worth a code ticket?

Wayfinder ticket: [Indel arm: is position-tolerant matching worth a code ticket?](https://github.com/dimadatascience/NanoGen/issues/91), map [Wayfinder: confirm genotype-v3 vs v3-ub on new diagnosis samples](https://github.com/dimadatascience/NanoGen/issues/72). Scored against §8 of [`confirmation-preregistration.md`](https://github.com/dimadatascience/NanoGen/blob/e6f629f/docs/research/confirmation-preregistration.md) (`e6f629f`), after the primary verdict was written ([`confirmation-verdict.md`](https://github.com/dimadatascience/NanoGen/blob/c637110/docs/research/confirmation-verdict.md), `c637110`). This arm cannot change that verdict.

**Patient data.** Aggregate numbers only. Targets are T1…Tn per patient (key in the HPC working folder's README).

## Answer

**Not worth a code ticket** under the pre-registered rule. Criterion 2 (no calibration harm) fails because the prototype fails the pooled calibration test (§5.1).

| §8 criterion | result | numbers |
|---|---|---|
| 1. Gain | **pass** | gated tumour MUT calls 1,210 → 1,922: **+712 (+58.8 %)** vs summed tie bands 195 and the 10 % bar |
| 2. No calibration harm | **fail** | prototype pooled O/E **1.66**, p **< 0.0001** (5 % of seeds pass); no new site exception at an indel target |
| 3. Precision | **pass** | FDR ≤ 0.2 at all 3 gated indel targets (max 0.14) |

**How to read the failure.** The criterion-2 failure is mostly **inherited, not caused**.
- The prototype's SNV rows are identical to v3's in all 8 patients, so its pool carries v3's own miss unchanged. v3 fails §5.1 by the same margin (O/E 1.66, p 0.0002; [#89](https://github.com/dimadatascience/NanoGen/issues/89)).
- The verdict kept v3 under rule 4 without a §5.1 pass. Once that happened, criterion 2 as written could not pass for any indel rule built on v3. The pre-registration did not foresee this case (it assumed that the rule §7 chose had passed). Reading criterion 2 as "no harm *relative to v3*" would be a post-hoc change, so it is **not** applied; the verdict stands as pre-registered.

**The arm is not a clean "harmless gain" either.** Relative to v3, the non-gating evidence shows real harm at two indel targets outside the gain set:
- **PALM08E T7:** negatives' ALT UMIs go from 0 to 10 (of 419), and tumour calls go from 0 to 45 at FDR 2.1. Under the prototype's own counts the target would be outlier-flagged (c_s 0.024 ≥ 0.5 × q_T 0.012).
- **PALM21E T7:** negatives' ALT UMIs go from 10 to 28, c_s from 0.025 to 0.077, and tumour FDR from 0.11 to 0.16.
- The 10-target indel sub-pool O/E rises from 1.35 to 1.50 (p 0.28 → 0.11). It still passes, but it is non-gating.

The gain is concentrated: **PALM11E T1** contributes 697 of the 712 extra calls (938 → 1,635). PALM08E T3 adds 15 and PALM14E T3 adds 0.

## Method

- **Call rule:** v3 (§7 rule 4): θ = min(1, P_s + c), m_min = 1, MUT iff m ≥ 1 and p < 0.1. `call` is imported from `bin/long_reads/genotyping.py` @ `61d8809`. The scoring library is the #89 code @ `23e407c` (`common89.py`, `score_lib.py`, and `holm`/`pooled` from `score89.py`), unmodified.
- **Contender:** the indel-arm run (`prototype/position-tolerant-indel`@`c77e147`, `--indel_window 12`, [#79](https://github.com/dimadatascience/NanoGen/issues/79)), with P_s, c_A, c_s and calls from its own counts. **Baseline:** the same rule on the v3 run.
- **Fixed sets:** the 57-target pool, the 20 seeds and the held-out split are those of §4–§5. The outlier flags and the no-signal list are the frozen #89 stage-1 tables (v3 run, §5.4, §6.2).
- **Code:** [`indel-arm/score91.py`](indel-arm/score91.py) and [`indel-arm/diag91.py`](indel-arm/diag91.py), committed at `e8f0afd` before the jobs ran. The jobs ran on `hpc_new` via sbatch + singularity in the pinned v3 image ([`run91.sbatch`](indel-arm/run91.sbatch)): job 9170839 (score) and job 9170840 (diagnostics).

**Checks (all pass, asserted in code):**
- The seeds, pool and tumour-evaluable set match the pre-registration.
- The pool's indel targets in the run outputs are exactly §8's 10. AML5 T5 is an indel but is outside the pool.
- Re-calling the prototype's counts with its own in-sample P_null reproduces its `genotype` column, with 0 differences in all 8 patients.
- The v3 arm reproduces #89's v3 pooled O, E, p and O/E exactly, for both the primary test and the fragility rerun.
- The v3 site exceptions are still exactly {AML5 T11}.

## Criterion 1: gain

**Gain set** (§8): PALM08E T3, PALM11E T1, PALM14E T3, sAML1 T11. **sAML1 T11 leaves the gate**: it is outlier-flagged (c_s 0.61, [#89](https://github.com/dimadatascience/NanoGen/issues/89)), and flagged targets are excluded from every tumour criterion (§5.4).

| target | status | O_T v3 | O_T prototype | FDR v3 | FDR prototype | tie band |
|---|---|---|---|---|---|---|
| PALM08E T3 | gated | 246 | 261 | 0.137 | 0.141 | 116 |
| PALM11E T1 | gated | 938 | 1,635 | 0.019 | 0.019 | 31 |
| PALM14E T3 | gated (thin null) | 26 | 26 | 0.089 | 0.090 | 48 |
| sAML1 T11 | flagged | 133 | 75 | 0.059 | 0.055 | 66 |
| **gated total** | | **1,210** | **1,922** | | | **195** |

**Tie bands:** band_s = max(⌈40 × rows_s / 34,979⌉, measured_s).
- measured_s = 0 in every patient. There are 0 SNV-row differences in WT, MUT or call across 367,705 SNV rows, so each band is its floor. The floors are PALM07E 38, PALM08E 116, PALM11E 31, PALM14E 48, PALM21E 137, AML4 25, AML5 51, sAML1 66.
- The band is summed over the patients with a gated target: 116 + 31 + 48 = 195.

**Sensitivity (non-gating):** including the flagged sAML1 T11 and its band gives 1,343 → 1,997 calls (+654) vs 261. That still passes, although tolerant matching **loses** 58 calls at sAML1 T11 itself.

**Phased locus:** PALM14E T3 is §8's named indel member of {T3–T5}. T4 and T5 are SNVs, identical across runs, so the §6.1 median-of-three rule gives the same zero gain.

## Criterion 2: calibration

Pooled Poisson-binomial (§5.1), seed medians, over the 37 unflagged pool targets (34 target-seed units per seed at the median):

| test | run | O | E | O/E (5–95 %) | p upper | p lower | pass | seeds passing |
|---|---|---|---|---|---|---|---|---|
| primary | v3 | 57.5 | 34.3 | 1.66 (1.22–2.20) | 0.0002 | 0.9999 | **fail** | 10 % |
| primary | prototype | 66.0 | 39.1 | 1.66 (1.23–2.10) | < 0.0001 | 1.0000 | **fail** | 5 % |
| fragility (n_A ≥ 11) | v3 | 55.0 | 34.3 | 1.61 | 0.0004 | | fail | 15 % |
| fragility (n_A ≥ 11) | prototype | 62.5 | 39.0 | 1.62 | 0.0002 | | fail | 15 % |
| indel sub-pool (reported) | v3 | 5.0 | 4.0 | 1.35 (0.48–2.48) | 0.28 | | (pass) | 65 % |
| indel sub-pool (reported) | prototype | 12.5 | 8.7 | 1.50 (0.89–2.12) | 0.11 | | (pass) | 50 % |

- The prototype adds ≈ 8.5 null calls in B against ≈ 4.8 more expected. All of the change is at indel targets, mainly PALM08E T7 (B calls 0 → 3.5, E_B 0.02 → 2.8) and PALM21E T7 (4 → 7, E_B 3.0 → 4.5).
- **Site screen (Holm α = 0.10):** the prototype's only site exception is AML5 T11 (an SNV), the same as v3's. **There is no new site exception at an indel target.** The closest is PALM11E T1 (B calls 0 → 0.5, O/E 8.9 on E_B 0.43; Holm-adjusted p 1.0).
- **Not fragile:** the primary test and the fragility rerun agree for both runs.

## Criterion 3: precision

FDR ≤ 0.2 at every gated clonal indel target: PALM08E T3 0.141, PALM11E T1 0.019, PALM14E T3 0.090. **Pass.**

## All pool indel targets (in-sample, both runs)

| target | flagged | ALT/total UMIs in negatives (v3 → proto) | q_T (v3 → proto) | tumour calls (v3 → proto) | FDR (v3 → proto) |
|---|---|---|---|---|---|
| PALM08E T3 | no | 10/94 → 10/90 | 0.25 → 0.28 | 246 → 261 | 0.14 → 0.14 |
| PALM08E T6 | no | 0/1 → 0/1 | 0.06 → 0.07 | 7 → 8 | 0.00 → 0.00 |
| PALM08E T7 | no (proto's own counts would flag it) | 0/466 → 10/419 | 0.00 → 0.01 | 0 → 45 | ∞ → **2.11** |
| PALM11E T1 | no | 0/695 → 1/688 | 0.06 → 0.12 | 938 → 1,635 | 0.02 → 0.02 |
| PALM14E T3 | no | 1/6 → 1/6 | 0.48 → 0.50 | 26 → 26 | 0.09 → 0.09 |
| PALM21E T1 | yes | 81/133 → 120/151 | 0.88 → 0.94 | 10 → 2 | 0.11 → 0.09 |
| PALM21E T2 | no | 0/226 → 0/218 | 0.00 → 0.01 | 1 → 2 | 0.03 → 0.01 |
| PALM21E T7 | no | 10/408 → 28/364 | 0.35 → 0.43 | 72 → 83 | 0.11 → **0.16** |
| AML5 T2 | no | 0/6 → 0/6 | 0 → 0 | 0 → 0 | – |
| sAML1 T11 | yes | 125/204 → 133/194 | 0.99 → 1.00 | 133 → 75 | 0.06 → 0.06 |

## Tolerant-match diagnostics (§10)

The diagnostics are summed over QC cells from the prototype's per-cell CONSENSUS tables. Only successful tasks are read, and each barcode once (each table also appears as a symlink staged into the aggregation task). Full table: [`out/diag_by_set.tsv`](indel-arm/out/diag_by_set.tsv).

| target | tol_reads (tumour / neg) | tol_unseen_reads | brc_mismatch_umis |
|---|---|---|---|
| PALM08E T3 | 105,874 / 531 | 0 | 0 |
| PALM08E T6 | 4,068 / 6 | 0 | 0 |
| PALM08E T7 | 581,025 / 33,571 | 0 | 0 |
| PALM11E T1 | 23,072 / 405 | 2 | 0 |
| PALM14E T3 | 23,244 / 39 | 0 | 0 |
| PALM21E T1 | 14,083 / 3,449 | 0 | 0 |
| PALM21E T2 | 7,119 / 37,037 | 1 | 5 |
| PALM21E T7 | 26,923 / 32,336 | 0 | 1 |
| AML5 T2 | 200 / 220 | 0 | 0 |
| sAML1 T11 | 75,270 / 19,251 | 17 | 0 |

- **The read classifier agrees with bam-readcount almost everywhere:** 6 mismatching UMIs in total, and 20 tolerant reads unseen at the anchor.
- **Tolerant reads are dominated by WT-majority molecules.** At PALM08E T7, 581 k tolerant tumour reads yield 85 MUT UMIs. At PALM21E T2 and T7 the negatives carry as many tolerant reads as the tumour cells. In other words, the ±12 bp / length ±1 window often matches non-somatic indel-like noise near those sites.

## Deviations

1. **Criterion 2 is applied literally,** although it cannot pass for any v3-based indel rule once v3 itself fails §5.1 (see Answer). The relative comparison is reported as non-gating.
2. **Outlier flags come from the v3 run** (§5.4) for both runs. The prototype's own counts would additionally flag PALM08E T7. That target is outside the gain set, so this is reported only.
3. **Duplicate tables:** every barcode's table is present twice (the task output plus a symlink staged into the aggregation task). The first successful copy is read, and the two copies are identical by construction.
4. **Denominators** are the cells in the outputs, as in #89 deviation 8 (the restart BAMs lack 10, 71 and 704 QC cells in AML4, AML5 and sAML1).

## Outputs

[`indel-arm/out/`](indel-arm/out/): `verdict.json`, `gain.tsv`, `pooled.tsv`, `pooled_by_seed.tsv`, `site_screen.tsv`, `targets.tsv`, `diag.tsv`, `diag_by_set.tsv`, `diag_meta.tsv`. HPC: `q91_indel/` (code, out, logs), with a "ticket 91" section in the README.
