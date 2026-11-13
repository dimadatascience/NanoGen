# Scoring: v3 vs v3-ub on the new diagnosis pool

Wayfinder ticket [Score v3 vs v3-ub against the confirmation pre-registration](https://github.com/dimadatascience/NanoGen/issues/89), map [Wayfinder: confirm genotype-v3 vs v3-ub on new diagnosis samples](https://github.com/dimadatascience/NanoGen/issues/72).

- **Pass bar:** [`confirmation-preregistration.md`](https://github.com/dimadatascience/NanoGen/blob/e6f629f/docs/research/confirmation-preregistration.md) at **`e6f629f`**, §4–§7 and §10. Deviations are listed [below](#deviations-and-interpretations).
- **Runs:** [Full runs: genotype-v3 and indel prototype on the new diagnosis samples](https://github.com/dimadatascience/NanoGen/issues/79), v3 @ `61d8809`, 8 patients.
- **Scoring code:** `docs/research/confirmation-score/` on this branch.
  - The no-signal list and the outlier flags were frozen at **`155b5cd`**, before either contender's tumour outcomes were computed (§6.2).
  - The scoring library is frozen at **`23e407c`**. `score_lib.py` is the P64 library, unmodified, and `call` is imported from `bin/long_reads/genotyping.py` @ `61d8809`.
  - Python is numpy 2.0.1, pandas 2.2.2 and scipy 1.14.0, from the pinned v3 image, on `hpc_new` via sbatch.
- **Patient data:** aggregate numbers only. Targets are T1…Tn, keyed in the HPC working folder's README.

The indel arm (§8) is not scored here. It runs after the verdict, in [Indel arm: is position-tolerant matching worth a code ticket?](https://github.com/dimadatascience/NanoGen/issues/91).

## Answer in brief

| criterion | v3 | v3-ub |
|---|---|---|
| §5.1 pooled calibration (median p ≥ 0.05 and O/E ≤ 1.5) | **fail**: O/E 1.66, p 0.0002, 10 % of seeds pass | **pass**: O/E 0.63, p 1.00, all seeds pass |
| §5.1 fragility rerun (n_A ≥ 11) | fail (O/E 1.61, p 0.0004) | pass (O/E 0.63) |
| §6.3 FDR ≤ 0.2 at every gated target | **pass** (max 0.17) | **pass** (max 0.09) |
| §6.5 tumour side evaluable | yes: 8 gated loci in 4 patients, 0 no-signal | same |
| §6.4 v3-ub tumour-call loss ≤ 5 % | n/a | **fail: 11.2 %** |
| **passes (§7)** | **no** (pooled calibration) | yes, but the loss cap blocks replacement |

**The verdict implied by §7 is rule 4: neither contender is confirmed, and v3 stays provisional.**
- v3 fails the pooled calibration, so rule 1 does not apply.
- v3-ub passes but loses 11.2 % of gated tumour calls, more than the 5 % cap, so rule 2 does not apply.

**Labels:**
- **fragile:** no. Both contenders keep their outcome under the n_A ≥ 11 rerun.
- **site exceptions:** AML5 T11 (v3 only); see §5.3.
- **weak-null pass:** yes, at PALM14E {T3–T5}; see §6.3.
- **blind spot (§1, verbatim):** "On this census, neither the pooled test nor the site screen can rule out a single-site miss the size of PALM09E T3's (O/E ≈ 1.7). A pooled pass therefore means 'no average miss ≥ 1.5× and no gross site miss on well-covered targets', **not** 'no T3-like site misses'."
- **scope line (§6.6):** "Per-target precision confirmed on 8 covered, signal-bearing loci in 4 patients. 24 under-covered and 0 no-signal targets are listed. Thin-null loci: PALM14E {T3–T5}. Weak-null pass: PALM14E {T3–T5}."

## Implementation checks (§3): pass

- **Re-calls reproduce the runs.** In all 8 patients, re-calling each v3 run's counts with its own `P_null`, using the imported `call`, reproduces the `genotype` column exactly. Recomputing c_s from the negatives matches the run to ≤ 1.1e-16.
- **Seeds and fixed sets are asserted in code.** The seed list is re-derived from `random.Random(20261112)`. The 57 pool targets are checked against `pool86.tsv` @ `546386e`, and the 17 tumour-evaluable targets against `leuk_covered ≥ 100` on the WES-clonal sets.
- **Cell sets** are QC cells only. Negatives (T/B/NK) are 3,351, 877, 293, 1,064, 4,431, 767, 518 and 1,073 (PALM07E, 08E, 11E, 14E, 21E, AML4, AML5, sAML1), with no overlap between negatives and Leukemic cells.

## Outlier flag (§5.4): 20 of 57 pool targets flagged

The flag is in-sample c_s ≥ 0.25 **or** c_s ≥ 0.5 × q_T, from the v3 run, with c_s = q_T = 0 unflagged. It removes **20 pool targets** from the pooled sums, the site screen and every tumour criterion.
- **13 flagged by c_s ≥ 0.25:** PALM07E T4; PALM14E T2 (c_s 0.94); PALM21E T1; AML4 T4; AML5 T1, T4, T10; sAML1 T7–T12.
- **7 flagged by the q_T clause alone:** PALM07E T2, PALM11E T3, PALM21E T5, AML4 T3, AML4 T9, AML5 T16 and sAML1 T2. At three of these (PALM07E T2, AML4 T3, AML4 T9), q_T = 0 with c_s > 0. The pre-registered rule flags this case literally; only c_s = q_T = 0 is exempt.

**Tumour-evaluable targets lost to the flag: 7 of 17.**
- PALM14E T2.
- AML5 T1.
- **All five of sAML1's** (T7, T9, T10, T11, T12; c_s 0.26–0.61).

sAML1 therefore keeps no gated locus. Its negatives read MUT at 26–61 % of UMIs at its tumour targets. That pattern fits negatives contaminated with tumour material, or a T/B/NK list that includes tumour cells, in this older-study sample. This report does not test either explanation.

The pooled calibration test runs on the remaining **37 targets in 8 patients**.

**Molecule purity** (descriptive) is not reported; see deviation 1.

## Calibration

### 5.1 Pooled calibration test

The test runs over every covered B cell at every included target-seed. Values are seed medians, with the 5–95 % range in parentheses.

| test | contender | O | E | O / E | p = P(S ≥ O) | P(S ≤ O) | target-seeds / seed | pass | seeds passing |
|---|---|---|---|---|---|---|---|---|---|
| primary | v3 | 57.5 (48.7–68.5) | 34.3 (29.9–40.8) | **1.66** (1.22–2.20) | **0.0002** (6.7e-09–0.082) | 1.00 | 34 | **fail** | 0.10 |
| primary | v3-ub | 38.5 (28.8–45) | 60.2 (54.5–66.7) | **0.63** (0.46–0.77) | 1.00 (0.97–1) | **0.0011** | 34 | **pass** | 1.00 |
| fragility (n_A ≥ 11) | v3 | 55 (44–67.6) | 34.3 (29.9–40.7) | 1.61 (1.20–2.13) | 0.0004 (3.6e-08–0.10) | 1.00 | 25 | fail | 0.15 |
| fragility (n_A ≥ 11) | v3-ub | 38.5 (28.8–45) | 60.2 (54.4–66.7) | 0.63 (0.46–0.77) | 1.00 (0.97–1) | 0.0011 | 25 | pass | 1.00 |

- **v3 over-calls in negatives** by about two-thirds on average. That exceeds the 1.5× bar, and only 2 of 20 seeds pass.
- **v3-ub over-corrects.** It makes 37 % fewer false calls than its own nominal expectation (lower-tail p = 0.001). The pre-registration does not gate on this; it prices it through the §6.4 loss cap.
- **Not fragile.** Dropping target-seeds with n_A < 11 leaves both outcomes unchanged.
- **Zero-information target-seeds** (n_A = 0, dropped for both contenders): **55 of 740 (7.4 %)**, under the 10 % warning line. By patient: PALM08E 17/120, PALM14E 12/120, AML4 11/80, sAML1 7/40, AML5 6/80, PALM21E 2/180, PALM07E 0/60 and PALM11E 0/60.

### 5.2 Per patient (non-gating)

Values are seed-median Σ O / Σ E. The exact 90 % Garwood interval on O is divided by E and reported as the seed medians of its bounds. LOPO is the pooled test with that patient left out (seed-median O/E and p).

| patient | v3 O / E | v3 O/E (90 % CI) | v3-ub O / E | v3-ub O/E (90 % CI) | LOPO v3 | LOPO v3-ub |
|---|---|---|---|---|---|---|
| PALM07E | 18 / 13.2 | 1.30 (0.82–1.99) | 13 / 17.3 | 0.76 (0.45–1.25) | 1.85, p 0.0002 | 0.55, p 1.00 |
| PALM08E | 3 / 2.19 | 1.43 (0.41–3.61) | 2.5 / 6.26 | 0.40 (0.09–1.12) | 1.68, p 0.0003 | 0.67, p 1.00 |
| PALM11E | 0 / 0.33 | 0 (0–9.0) | 0 / 5.11 | 0 (0–0.59) | 1.68, p 0.0002 | 0.69, p 0.99 |
| PALM14E | 13 / 3.90 | **2.90** (1.72–4.58) | 7 / 3.21 | 1.96 (0.95–3.61) | 1.52, p 0.003 | 0.56, p 1.00 |
| PALM21E | 11.5 / 8.25 | 1.37 (0.84–2.22) | 9 / 22.6 | 0.41 (0.22–0.71) | 1.79, p 0.0002 | 0.79, p 0.92 |
| AML4 | 7 / 4.04 | 1.62 (0.73–3.08) | 5 / 5.52 | 0.99 (0.41–2.03) | 1.65, p 0.0002 | 0.61, p 1.00 |
| AML5 | 3 / 0.66 | **5.06** (1.38–13.0) | 2 / 1.06 | 2.06 (0.37–6.69) | 1.59, p 0.001 | 0.59, p 1.00 |
| sAML1 | 1.5 / 0.74 | 1.68 (0.20–7.97) | 0 / 0.19 | 0 (0–15.7) | 1.64, p 0.0004 | 0.64, p 1.00 |

- **v3's miss is not driven by one patient.** Its point O/E exceeds 1 in 7 of 8 patients, and every leave-one-patient-out pool still fails (O/E 1.52–1.85, p ≤ 0.003).
- **The 90 % CI excludes 1 only at PALM14E and AML5.** The excess concentrates where the A-half estimate c_A is small but B still reads MUT: PALM14E T6, AML5 T3 and AML5 T11.

### 5.3 Site screen (flags, never vetoes)

The screen runs on the 37 unflagged pool targets. Each target's p is the per-seed exact Poisson tail, taken as the median across seeds; Holm is applied at α = 0.10 separately for each contender. Every target had at least one included seed, so none is "no negative data".

**Site exception: AML5 T11, for v3 only.**
- **v3:** O_B 1 against E_B 0.0009 (median p 0.001, Holm-adjusted 0.037).
- **v3-ub:** O_B 0 against E_B 0, unflagged.
- **Cause:** the target has 11 covered negatives, with one MUT UMI in-sample (c_s 0.077). In most seeds A holds no MUT, so v3's held-out θ falls to P_s and a single ALT UMI in B is called MUT.
- **Effect:** the target stays in Σ O / Σ E and does not change the decision rule. Per §5.3 it opens a germline/artefact follow-up. That is graduated by [Verdict: does genotype-v3 hold, or does v3-ub replace it?](https://github.com/dimadatascience/NanoGen/issues/90).

**Next-smallest v3 site p** (none survives Holm):
- PALM14E T6: 7 vs 2.2, p 0.009, Holm 0.32.
- AML5 T3: 3 vs 0.63, p 0.035.
- No v3-ub target is below 0.03 before correction.

### Per-target calibration table (§10)

Values are seed medians. Status gives the outcome (tested, site exception or flagged); flagged rows are reported but excluded from every sum.

| patient | target | neg covered | B covered | A UMIs | c_s | CPU(x_s, n_s) | c_A | q_T | seeds incl. | v3 O_B / E_B | v3 site p | v3-ub O_B / E_B | v3-ub site p | status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PALM07E | T1 | 279 | 140 | 200 | 0.027 | 0.044 | 0.024 | 0.46 | 20 | 6 / 4.3 | 0.24 | 3.5 / 6.9 | 0.93 | tested |
| PALM07E | T2 | 13 | 6 | 8 | 0.067 | 0.28 | 0.095 | 0 | 20 | 0 / 0.0008 | 1 | 0 / 0 | 1 | flagged |
| PALM07E | T3 | 74 | 37.5 | 56 | 0.043 | 0.089 | 0.042 | 0.36 | 20 | 1.5 / 1.6 | 0.71 | 1 / 0.39 | 0.31 | tested |
| PALM07E | T4 | 12 | 5 | 9 | 0.47 | 0.7 | 0.5 | 0.4 | 20 | 0 / 0 | 1 | 0 / 0 | 1 | flagged |
| PALM07E | T5 | 336 | 168 | 323 | 0.049 | 0.065 | 0.049 | 0.21 | 20 | 11 / 7.3 | 0.16 | 8.5 / 9.4 | 0.69 | tested |
| PALM08E | T1 | 6 | 3 | 3 | 0 | 0.39 | 0 | 0.54 | 20 | 0 / 0.0003 | 1 | 0 / 0 | 1 | tested |
| PALM08E | T3 | 66 | 33 | 45.5 | 0.11 | 0.17 | 0.11 | 0.25 | 20 | 0 / 0.22 | 1 | 0 / 0.3 | 1 | tested |
| PALM08E | T4 | 117 | 58 | 199 | 0.02 | 0.037 | 0.016 | 0.17 | 20 | 2.5 / 1.8 | 0.37 | 2.5 / 2.3 | 0.54 | tested |
| PALM08E | T5 | 43 | 21 | 45.5 | 0 | 0.032 | 0 | 0.0006 | 20 | 0 / 0.0047 | 1 | 0 / 1.1 | 1 | tested |
| PALM08E | T6 | 1 | 1 | 0 | 0 | 0.95 | 0 | 0.063 | 3 | 0 / 0 | 1 | 0 / 0 | 1 | tested |
| PALM08E | T7 | 166 | 83 | 231 | 0 | 0.0064 | 0 | 0 | 20 | 0 / 0.024 | 1 | 0 / 2.5 | 1 | tested |
| PALM11E | T1 | 225 | 115 | 350 | 0 | 0.0043 | 0 | 0.058 | 20 | 0 / 0.32 | 1 | 0 / 3.1 | 1 | tested |
| PALM11E | T3 | 12 | 5.5 | 6.5 | 0.17 | 0.44 | 0.13 | 0.27 | 20 | 0 / 0 | 1 | 0 / 0 | 1 | flagged |
| PALM11E | T4 | 21 | 10.5 | 34 | 0 | 0.041 | 0 | 0 | 20 | 0 / 0.0037 | 1 | 0 / 0.44 | 1 | tested |
| PALM11E | T5 | 100 | 49 | 105 | 0 | 0.015 | 0 | 0 | 20 | 0 / 0.01 | 1 | 0 / 1.5 | 1 | tested |
| PALM14E | T1 | 14 | 6.5 | 14 | 0 | 0.12 | 0 | 0.41 | 20 | 0 / 0.0009 | 1 | 0 / 0.029 | 1 | tested |
| PALM14E | T2 | 15 | 7 | 9 | 0.94 | 1 | 0.92 | 1 | 20 | 0 / 0 | 1 | 0 / 0 | 1 | flagged |
| PALM14E | T3 | 4 | 2 | 3 | 0.17 | 0.58 | 0.18 | 0.48 | 16 | 0 / 0.0002 | 1 | 0 / 0 | 1 | tested |
| PALM14E | T4 | 4 | 2 | 3 | 0.17 | 0.58 | 0.18 | 0.38 | 16 | 0 / 0.0002 | 1 | 0 / 0 | 1 | tested |
| PALM14E | T5 | 4 | 2 | 3 | 0.17 | 0.58 | 0.18 | 0.43 | 16 | 0 / 0.0002 | 1 | 0 / 0 | 1 | tested |
| PALM14E | T6 | 116 | 58 | 149 | 0.055 | 0.082 | 0.052 | 0.51 | 20 | 7 / 2.2 | 0.009 | 7 / 2.6 | 0.033 | tested |
| PALM14E | T7 | 82 | 40 | 93 | 0.098 | 0.14 | 0.1 | 0.51 | 20 | 2 / 1.5 | 0.66 | 0 / 0.51 | 1 | tested |
| PALM21E | T1 | 76 | 36 | 68.5 | 0.61 | 0.68 | 0.64 | 0.88 | 20 | 0 / 0.068 | 1 | 0 / 0 | 1 | flagged |
| PALM21E | T2 | 91 | 44.5 | 122 | 0 | 0.013 | 0 | 0.0039 | 20 | 0 / 0.01 | 1 | 0 / 1.3 | 1 | tested |
| PALM21E | T3 | 3 | 1.5 | 1.5 | 0 | 0.63 | 0 | 0.083 | 18 | 0 / 0.0001 | 1 | 0 / 0 | 1 | tested |
| PALM21E | T4 | 34 | 18.5 | 15.5 | 0 | 0.084 | 0 | 0.37 | 20 | 0 / 0.14 | 1 | 0 / 0 | 1 | tested |
| PALM21E | T5 | 18 | 9 | 14 | 0.18 | 0.34 | 0.17 | 0.33 | 20 | 0 / 0.11 | 1 | 0 / 0.042 | 1 | flagged |
| PALM21E | T7 | 185 | 93 | 204 | 0.025 | 0.041 | 0.024 | 0.35 | 20 | 4 / 3 | 0.34 | 4 / 4 | 0.52 | tested |
| PALM21E | T8 | 282 | 138 | 184 | 0.017 | 0.032 | 0.015 | 0.081 | 20 | 3 / 2.4 | 0.46 | 3 / 5.5 | 0.93 | tested |
| PALM21E | T9 | 122 | 61.5 | 62.5 | 0.1 | 0.16 | 0.11 | 0.35 | 20 | 0 / 0.043 | 1 | 0 / 0.065 | 1 | tested |
| PALM21E | T10 | 228 | 118 | 116 | 0.0083 | 0.026 | 0.0086 | 0.12 | 20 | 1 / 1.1 | 0.66 | 1 / 4.9 | 0.99 | tested |
| PALM21E | T11 | 493 | 244 | 276 | 0.0036 | 0.011 | 0.0037 | 0.053 | 20 | 1 / 1 | 0.65 | 1 / 4.8 | 0.99 | tested |
| PALM21E | T12 | 137 | 70 | 163 | 0 | 0.0092 | 0 | 0 | 20 | 0 / 0.016 | 1 | 0 / 1.8 | 1 | tested |
| AML4 | T1 | 107 | 54 | 107 | 0.18 | 0.23 | 0.17 | 0.86 | 20 | 2 / 0.8 | 0.39 | 1 / 0.88 | 0.51 | tested |
| AML4 | T2 | 169 | 86 | 116 | 0.037 | 0.063 | 0.031 | 0.73 | 20 | 4.5 / 3.2 | 0.32 | 4 / 4.6 | 0.75 | tested |
| AML4 | T3 | 5 | 2 | 4.5 | 0.12 | 0.47 | 0.14 | 0 | 18 | 0 / 0.0004 | 1 | 0 / 0 | 1 | flagged |
| AML4 | T4 | 23 | 12 | 14.5 | 0.48 | 0.64 | 0.45 | 0.7 | 20 | 0 / 0.035 | 1 | 0 / 0 | 1 | flagged |
| AML4 | T6 | 8 | 4 | 6 | 0.091 | 0.36 | 0.056 | 1 | 20 | 0.5 / 0.0006 | 0.5 | 0 / 0 | 1 | tested |
| AML4 | T9 | 13 | 7 | 8.5 | 0.056 | 0.24 | 0.077 | 0 | 20 | 0 / 0.044 | 1 | 0 / 0.019 | 1 | flagged |
| AML4 | T10 | 1 | 1 | 0 | 0 | 0.78 | 0 | 0 | 9 | 0 / 0 | 1 | 0 / 0 | 1 | tested |
| AML5 | T1 | 56 | 28.5 | 38.5 | 0.46 | 0.56 | 0.49 | 1 | 20 | 1 / 0.072 | 0.036 | 1 / 0.063 | 0.067 | flagged |
| AML5 | T2 | 2 | 1 | 5 | 0 | 0.39 | 0 | 0 | 16 | 0 / 0.0001 | 1 | 0 / 0 | 1 | tested |
| AML5 | T3 | 147 | 74 | 99 | 0.19 | 0.24 | 0.18 | 0.55 | 20 | 3 / 0.63 | 0.035 | 2 / 1.1 | 0.25 | tested |
| AML5 | T4 | 11 | 5.5 | 16 | 0.74 | 0.86 | 0.74 | 0.79 | 20 | 0 / 0 | 1 | 0 / 0 | 1 | flagged |
| AML5 | T7 | 2 | 1 | 2 | 0 | 0.63 | 0 | 0 | 18 | 0 / 0.0001 | 1 | 0 / 0 | 1 | tested |
| AML5 | T10 | 1 | 0.5 | 0.5 | 1 | 1 | 1 | 0.6 | 10 | 0 / 0 | 1 | 0 / 0 | 1 | flagged |
| AML5 | T11 | 11 | 6 | 6 | 0.077 | 0.32 | 0 | 0.53 | 20 | 1 / 0.0009 | 0.001 | 0 / 0 | 1 | **site exception (v3)** |
| AML5 | T16 | 57 | 29.5 | 64 | 0.1 | 0.16 | 0.11 | 0.17 | 20 | 1 / 0.5 | 0.38 | 1 / 0.5 | 0.45 | flagged |
| sAML1 | T1 | 38 | 19 | 26.5 | 0.054 | 0.13 | 0.051 | 0.22 | 20 | 1.5 / 0.74 | 0.39 | 0 / 0.19 | 1 | tested |
| sAML1 | T2 | 6 | 3 | 3 | 0.17 | 0.58 | 0.23 | 0.22 | 20 | 0 / 0 | 1 | 0 / 0 | 1 | flagged |
| sAML1 | T4 | 2 | 1 | 1 | 0 | 0.78 | 0 | 0.53 | 13 | 0 / 0 | 1 | 0 / 0 | 1 | tested |
| sAML1 | T7 | 45 | 21.5 | 39 | 0.49 | 0.59 | 0.47 | 0.4 | 20 | 1 / 0.083 | 0.085 | 1 / 0.05 | 0.044 | flagged |
| sAML1 | T8 | 77 | 40 | 69.5 | 0.27 | 0.34 | 0.24 | 0.15 | 20 | 2 / 0.47 | 0.21 | 1 / 0.37 | 0.31 | flagged |
| sAML1 | T9 | 241 | 123 | 188 | 0.26 | 0.3 | 0.25 | 0.52 | 20 | 7 / 1.7 | 0.0025 | 6.5 / 1.9 | 0.012 | flagged |
| sAML1 | T10 | 51 | 24.5 | 44.5 | 0.54 | 0.63 | 0.55 | 0.69 | 20 | 1 / 0.08 | 0.034 | 1 / 0.031 | 0.077 | flagged |
| sAML1 | T11 | 103 | 51 | 104 | 0.61 | 0.67 | 0.62 | 0.99 | 20 | 2 / 0.24 | 0.014 | 0 / 0.071 | 1 | flagged |
| sAML1 | T12 | 24 | 12 | 26.5 | 0.37 | 0.49 | 0.32 | 0.44 | 20 | 1.5 / 0.15 | 0.036 | 0.5 / 0.07 | 0.56 | flagged |

`out/site_screen.tsv` and `out/targets.tsv` hold the full set: covered Ambiguous and Leukemic cells; ranges for c_A, CPU(x_A, n_A), O_B and E_B; lower-tail p; the false-MUT rate in B with Wilson 95 % CIs; and the Holm-adjusted p.

## Tumour side

### 6.2 No-signal list (frozen at `155b5cd`)

Every tumour-evaluable, unflagged target carries signal: **the no-signal list is empty.**
- The 10 such targets give π̂ 0.20–1.00 and f̂ ≥ ρ₀ + 0.05.
- The closest case is PALM11E T1: π̂ 1.00 and f̂ 0.058 against ρ₀ + 0.05 = 0.051.

### 6.1 / 6.5 Evaluability

There are **15 loci** in all.

| status | count | loci |
|---|---|---|
| gated | 8 | PALM08E T3, T4; PALM11E T1; PALM14E T1, {T3–T5}, T6, T7; AML5 T3 |
| flagged | 7 | PALM14E T2; AML5 T1; sAML1 T7, T9, T10, T11, T12 |
| no-signal | 0 | – |

- **Patients with a gated locus: 4** (PALM08E, PALM11E, PALM14E, AML5). That is ≥ 3, and 0 no-signal loci is under 3, so **the tumour side is evaluable** (§6.5).
- **Under-covered clonal targets (24),** listed and never counted:
  - AML4 T1–T4, T6–T10;
  - AML5 T2, T5–T7, T9, T10, T12;
  - PALM08E T1, T2;
  - sAML1 T2–T5, T13, T14.
- **Coverage drop below 100** (§6.1, kept in). In the v3 run, PALM14E T1 has 99 covered tumour cells and T3/T4/T5 have 97/95/94. The pre-registered screen (`pool86.tsv`) put them at ≥ 100.

### 6.3 Precision gate, 6.4 loss and efficiency

These are in-sample, over covered tumour cells. FDR = E_T / O_T, with the estimated true null ρ = min(1, P_s + c_s). Efficiency (P64 3b) uses the shared v3-count mixture fit. "Undet." is the undetermined fraction of covered tumour cells.

| patient | target | status | tumour covered | n_s | π̂ / f̂ | v3 O_T | v3 FDR | v3-ub O_T | v3-ub FDR | loss | efficiency v3 / ub | undet. v3 / ub |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PALM08E | T3 | gated | 2089 | 94 | 0.20 / 0.85 | 246 | 0.14 | 228 | 0.056 | 0.073 | 0.51 / 0.51 | 0.61 / 0.61 |
| PALM08E | T4 | gated | 1910 | 391 | 0.21 / 0.87 | 426 | 0.17 | 413 | 0.073 | 0.031 | 0.94 / 0.99 | 0.31 / 0.32 |
| PALM11E | T1 | gated | 3970 | 695 | 1.00 / 0.058 | 938 | 0.019 | 931 | 0.018 | 0.0075 | 0.24 / 0.23 | 0.16 / 0.16 |
| PALM14E | T1 | gated | 99 | 23 | 0.36 / 0.93 | 34 | 0.0006 | 7 | 1e-08 | 0.79 | 0.96 / 0.20 | 0.61 / 0.88 |
| PALM14E | T3 | gated (thin null) | 97 | 6 | 0.41 / 1.00 | 26 | 0.089 | 9 | 9e-05 | 0.65 | 0.63 / 0.23 | 0.54 / 0.70 |
| PALM14E | T4 | gated (thin null) | 95 | 6 | 0.32 / 0.98 | 17 | 0.081 | 4 | 0.0002 | 0.76 | 0.52 / 0.13 | 0.62 / 0.76 |
| PALM14E | T5 | gated (thin null) | 94 | 6 | 0.39 / 1.00 | 18 | 0.085 | 4 | 0.0001 | 0.78 | 0.47 / 0.11 | 0.63 / 0.78 |
| PALM14E | T6 | gated | 825 | 292 | 0.52 / 0.95 | 438 | 0.063 | 437 | 0.057 | 0.002 | 0.99 / 0.99 | 0.22 / 0.23 |
| PALM14E | T7 | gated | 636 | 174 | 0.43 / 0.95 | 282 | 0.16 | 118 | 0.027 | 0.58 | 0.95 / 0.43 | 0.42 / 0.67 |
| AML5 | T3 | gated | 881 | 200 | 0.94 / 0.57 | 152 | 0.12 | 111 | 0.091 | 0.27 | 0.18 / 0.13 | 0.77 / 0.81 |
| PALM14E | T2 | flagged | 117 | 17 | – | 0 | ∞ | 0 | ∞ | – | – | 0.04 / 0 |
| AML5 | T1 | flagged | 828 | 82 | – | 192 | 0.070 | 107 | 0.028 | 0.44 | – | 0.77 / 0.87 |
| sAML1 | T7 | flagged | 381 | 75 | – | 5 | 0.22 | 1 | 0.12 | 0.80 | – | 0.85 / 0.86 |
| sAML1 | T9 | flagged | 738 | 380 | – | 69 | 0.19 | 68 | 0.19 | 0.014 | – | 0.84 / 0.84 |
| sAML1 | T10 | flagged | 603 | 85 | – | 59 | 0.11 | 21 | 0.023 | 0.64 | – | 0.78 / 0.84 |
| sAML1 | T11 | flagged | 1084 | 204 | – | 133 | 0.059 | 97 | 0.032 | 0.27 | – | 0.88 / 0.91 |
| sAML1 | T12 | flagged | 154 | 52 | – | 12 | 0.10 | 10 | 0.034 | 0.17 | – | 0.75 / 0.76 |

**Precision gate (§6.3): both contenders pass at every gated target.**
- v3's FDR ranges from 0.0006 to 0.17, with the highest values at PALM08E T4 (0.17) and PALM14E T7 (0.16).
- v3-ub's maximum FDR is 0.091.
- Flagged rows are shown for completeness and are not gated.

**Weak-null check (§6.3) at PALM14E {T3–T5}: weak-null pass.**
- v3 passes under the point-estimate null (FDR 0.081–0.089).
- With the null raised to CPU(x_s, n_s) = 0.58 (x_s = 1 of n_s = 6 UMIs), v3's FDR is **1.12, 1.16 and 1.06** at T3, T4 and T5.
- So v3 passes there only under the point estimate.

**v3-ub tumour-call loss (§6.4): 11.2 %, fail (cap 5 %).**
- Σ O_T over the 8 gated loci is 2,534 for v3 against 2,249 for v3-ub. The phased locus enters once, as the median of T3–T5: 18 vs 4.
- **By patient:**
  - PALM08E 4.6 %;
  - PALM11E 0.7 %;
  - **PALM14E 26.7 %;**
  - **AML5 26.9 %.**
- **Where the loss sits.** It concentrates at targets with few negative UMIs, where CPU(x_s, n_s) far exceeds c_s:
  - PALM14E T1: n_s 23, so CPU 0.12 at c_s = 0, and a 79 % loss;
  - PALM14E {T3–T5}: 65–78 %;
  - PALM14E T7: 58 %;
  - AML5 T3: 27 %.
- At well-sampled nulls the loss is small: PALM14E T6 0.2 %, PALM11E T1 0.7 %, PALM08E T4 3 %.
- **Deployment cost outside the gated set** (not gating): over the 30 other unflagged targets with covered tumour cells, v3-ub makes 601 tumour MUT calls against v3's 775 (−22 %). It silences 7 targets completely.

## Deviations and interpretations

1. **Molecule purity not reported** (§5.4, descriptive). The run outputs hold cell-level WT/MUT/OTH counts only, with no per-UMI read composition. Reporting it would need a new counting pass, outside this ticket. It is descriptive and changes nothing.
2. **CPU(x, 0).** The P64 library's `cp_upper` returns 0 at n = 0. Scoring uses the pre-registration's own definition (1 at n = 0 or x = n), in `common89.cpu`. This is moot: n_A = 0 target-seeds are dropped for both contenders, and n_s ≥ 1 at every target.
3. **Per-patient Garwood interval.** §5.2 does not say which seed carries the interval. It is computed per seed on (Σ O, Σ E), and the seed medians of its bounds are reported.
4. **Flagged targets and the site screen.** Flagged pool targets are excluded from the pooled sums (§5), so they are not among the "tested targets" of §5.3 and do not enter the Holm family (37 targets). Their site p is still listed.
5. **Literal q_T clause.** The flag is applied as written: c_s ≥ 0.5 q_T, with only c_s = q_T = 0 exempt. It therefore flags PALM07E T2, AML4 T3 and AML4 T9, where q_T = 0 < c_s. All three are calibration-only targets with ≤ 6 covered tumour cells and ≤ 0.044 expected null calls in B per seed. Unflagging them would add < 0.1 to v3's pooled E, which cannot move the verdict.
6. **Weak-null computation.** v3's calls are kept, with θ = P_s + c_s. Only the estimated true null in E_T is raised to min(1, P_s + CPU(x_s, n_s)).
7. **Phased-locus status.** The locus would count as no-signal only if all three members were. This is moot, because none is.
8. **Missing QC cells** (from the full-runs ticket). The restart BAMs lack 10, 71 and 704 QC barcodes (AML4, AML5, sAML1), so denominators are the cells in the outputs.
9. **Indel-arm quantities** (§10: tol_reads, tie bands) are deferred to the indel-arm ticket, which §8 schedules after the verdict.

## Reproduction

On `hpc_new`, the HPC folder is `q89_score/` (code/, out/, logs/; jobs 9170744 stage 1 and 9170747 stage 2, `short` partition, ≤ 1 CPU-h each):

```
sbatch --chdir=$PWD -o logs/stage1_%j.log code/run89.sbatch stage1   # flags + no-signal (commit before stage 2)
sbatch --chdir=$PWD -o logs/stage2_%j.log code/run89.sbatch stage2   # everything else
python3 tables89.py out                                               # tables in this report
```
