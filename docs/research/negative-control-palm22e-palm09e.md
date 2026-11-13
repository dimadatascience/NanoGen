# Scoring: T/B/NK negative controls vs blast calls (PALM22E, PALM09E)

Wayfinder ticket *Score the contenders against the pass bar* (#67), map *T/B negative controls vs blast calls* (#60). The pass bar, split, seeds, evaluability rules and computations are those of `docs/research/negative-control-preregistration.md` at **`76c3f99`**. Deviations are listed under [Deviations](#deviations-from-the-pre-registration). The scoring sections say only pass / fail / not evaluable. The decision on which contender NanoGen uses was made in *Verdict: which genotyping version to use for blast calls* (#69) and is recorded under [Verdict](#verdict).

**Patient data:** aggregate numbers only, targets as T1…Tn (key in the HPC working folder's README).

## Answer in brief

| contender | PALM22E | PALM09E |
|---|---|---|
| `main` (v1) | **fail**: 1 at T4, T5 (O/E 1,210 and 249); 2 at T4 (FDR 0.98) | **fail**: 1 at T1, T3; 2 at T1 (FDR 0.89) |
| v2 | **fail**: 1 at T4, T5 (O/E 29 and 48) | **fail**: 1 at T1, T3; 2 at T1 (FDR 0.25) |
| v2 + min2 | **pass** | **fail**: 1 at T3; 2 at T1 (O_T = 0, FDR ∞) |
| v3 | **pass** | **fail**: 1 at T3 (median O/E 1.70, p 0.048); 2 at T1 (FDR 0.25); 3a at T1 (tie, +34) |

- **PALM22E: v3 and v2 + min2 both pass.** v3 keeps about twice the tumour MUT calls at every clonal target: +1,457 at T5, +1,297 at T4, +763 at T3. That is far beyond the 40-call noise floor, so 3a passes and no replicate run is needed.
- **v3 is the only rule that calibrates without discarding single-UMI evidence.** Its median held-out O/E is 0.91 at T4 and 0.97 at T5. v2 + min2 passes calibration by making almost no calls (O_B = 0).
- **`main` and v2 are miscalibrated by 1–3 orders of magnitude in negatives.** `main`'s nominal null is the 1e-4 floor in most cells (median θᵢ = 1e-4), and at the two insertion targets of the same mutation type (PALM22E T4, PALM09E T1) `main`'s own counts read MUT in 49 % and 20 % of negative-cell UMIs.
- **PALM09E fails every contender.** Two causes, both already known:
  - **T1** (the sole clonal target) has almost no MUT signal in the v3 run's counts: q_T = 0.42 %, 34 calls in 1,428 cells. No call rule over these counts passes 2 or 3a. This is the counts issue that *PALM09E indel MUT loss in v3 consensus counts* (#70) investigates.
  - **T3**: v3 fails calibration on both conditions (median p 0.048 < 0.05, median O/E 1.70 > 1.5; 5 of 20 seeds pass). The cause is unresolved (see [Verdict](#verdict)). The exploratory upper-bound variant v3-ub passed it in the sweep (#68), but it was chosen on this target, so that result is not confirmatory.
- **Ties:** PALM09E fails every contender, so no tie decides anything there, and the PALM22E differences are all ≥ 763 calls. **No replicate runs** under the pre-registered rule.

## Inputs

- **Runs** from *Full runs: genotype-v3 and main on PALM22E and PALM09E* (#66):
  - genotype-v3 @ `61d8809`, `--negative_cells` = T/B/NK, `min_mut_umis = 1`, image `ghcr.io/dimadatascience/nanogen@sha256:764368d6…`;
  - `main` @ `9fe147a` through the tagged-BAM wrapper, `baseline_error = 1e-4`.
- **Cells:** QC cells only, matched on the bare 16-mer.
  - Negatives = T/B/NK list.
  - Tumour = Leukemic group; Ambiguous is reported only (from *Tumour compartment definition*, #61).
  - The run outputs contain 7,172 of 7,603 PALM22E QC cells and 3,159 of 4,966 PALM09E QC cells. These are the denominators.
- **Clonal targets** from *Clonal targets in PALM22E and PALM09E* (#63): PALM22E T5, T3, T4; PALM09E T1.
- **Software:** scoring in the pinned v3 image (numpy 2.0.1), Slurm job 9002995 on `short` (1 CPU, 16 s, 0.24 GB).
  - `call` is imported from `bin/long_reads/genotyping.py` @ `61d8809` (md5 `b3852efa…`) for every contender, including `main`.
  - `main`'s own `genotyping.py` @ `9fe147a` has md5 `cfd82f55…`.

### Implementation checks (all passed; `negative-control/score-out/checks_*.json`)

| check | PALM22E | PALM09E |
|---|---|---|
| v3: re-calling every row with the run's P_null reproduces `genotype` | ✓ | ✓ |
| v3: c_s recomputed over QC negatives vs the run | 3e-17 | 8e-17 |
| v3: re-calling QC rows with the recomputed c_s reproduces `genotype` | ✓ | ✓ |
| `main`: p recomputed from (m, k, θᵢ) by the imported `call`, max abs diff (m ≥ 1) | 1e-16 | 1e-16 |
| `main`: p is NaN wherever m = 0 | ✓ | ✓ |
| `main`: imported `call` with θᵢ reproduces `main`'s `genotype` exactly | ✓ | ✓ |
| `main`: one row per cell × target | ✓ | ✓ |

A synthetic dataset (3,000 cells, two targets, `main`'s output produced by `main`'s own `genotyping.py`) ran end to end first and passed the same checks.

## Evaluability and flags

Evaluable = ≥ 100 covered B cells and ≥ 200 WT + MUT UMIs in A (v3 counts), in ≥ 10 of 20 seeds.

| sample | target | clonal | covered Neg / Leuk (v3) | eval. seeds | c_s | c_s^main | q_T | flag |
|---|---|---|---|---|---|---|---|---|
| PALM22E | T1, T2, T6 | – | ≤ 10 / ≤ 42 | 0 | 0 | ≤ 0.009 | – | – |
| PALM22E | T3 | yes | 301 / 1,576 | 6 | 0.018 | 0.053 | 0.52 | – |
| PALM22E | T4 | yes | 2,029 / 4,237 | 20 | 0.006 | **0.489** | 0.17 | – |
| PALM22E | T5 | yes | 1,209 / 3,610 | 20 | 0.005 | 0.026 | 0.49 | – |
| PALM22E | T7 | – | 183 / 722 | 0 | 0.026 | 0.020 | 0.51 | – |
| PALM09E | T1 | yes | 956 / 1,428 | 20 | 0.0009 | **0.195** | 0.0042 | – |
| PALM09E | T2 | – | 678 / 785 | 20 | 0.0034 | 0.004 | 0.0034 | **c_s ≥ 0.5 q_T** |
| PALM09E | T3 | – | 312 / 220 | 20 | 0.040 | 0.024 | 0.13 | – |
| PALM09E | T4 | – | 147 / 139 | 0 | 0.007 | 0.020 | 0.09 | – |

- Criterion 1 is therefore decided on PALM22E T4, T5 and PALM09E T1, T3. PALM22E T3 is **not evaluable** for 1, but as a clonal target it still counts for 2 and 3a.
- PALM09E T2 is flagged: negatives and tumour read MUT at the same rate. It is scored below and excluded from the verdict.

## Criterion 1: held-out calibration (seed medians over evaluable seeds)

Each cell reads median O_B / E_B, then median p_cal, then the outcome, with the fraction of seeds passing in brackets. Pass = p_cal ≥ 0.05 and O/E ≤ 1.5.

| sample | target | covered B cells | `main` | v2 | v2 + min2 | v3 |
|---|---|---|---|---|---|---|
| PALM22E | T4 | ~1,010 | 1,210, 1e-14 ✗ (0) | 28.8, 0 ✗ (0) | 0 obs ✓ (0.55) | **0.91, 0.71 ✓ (0.85)** |
| PALM22E | T5 | ~606 | 249, 5e-15 ✗ (0) | 47.9, 1e-6 ✗ (0) | 0 obs ✓ (1.0) | **0.97, 0.59 ✓ (0.70)** |
| PALM09E | T1 | ~478 | 1,010, 8e-15 ✗ (0) | 5.85, 0.12 ✗ (0.2) | 0 obs ✓ (1.0) | **0.98, 0.64 ✓ (0.55)** |
| PALM09E | T3 | ~154 | 307, 2e-15 ✗ (0) | 402, 3e-16 ✗ (0) | 2.4e4, 4e-5 ✗ (0.35) | **1.70, 0.048 ✗ (0.25)** |
| PALM09E | T2 (flagged) | ~335 | 36.9 ✗ | 34.8 ✗ | ✓ | 1.03, 0.54 ✓ |

"0 obs" means O_B = 0: no MUT call in B, so the ratio is 0 (E_B > 0) and p = 1. The one-sided test passes this case by construction (see deviation 8).

- **5–95 % ranges** are wide where E_B is small. v3's O/E runs 0.73–1.63 at PALM22E T4, 0.48–3.4 at PALM22E T5 and 0–11.9 at PALM09E T1. Every column is in `score_<S>.tsv`.
- **Observed false-MUT rate in negatives, in-sample:**
  - v3: 2.7 % (PALM22E T4), 0.7 % (T5), 0.3 % (PALM09E T1), 7.4 % (PALM09E T3).
  - `main`: 95 %, 4.7 %, 66 %, 29 %.
- These match the pre-registration's expectation for v3: α = 0.1 per cell, so ~3–9 % false MUT per covered WT cell at c_s 1–5 %.

## Criterion 2 and sensitivity: tumour compartment, in-sample

ρ is the estimated true null:

- P_s + c_s for the v3-count contenders;
- θᵢ + c_s^main per cell for `main`.

FDR = E_T / O_T; pass ≤ 0.2. Efficiency = (O_T − (1 − π̂) E_T) / (π̂ n_cov,T), reported against 0.7 (3b, not gating).

| sample | target | metric | `main` | v2 | v2 + min2 | v3 |
|---|---|---|---|---|---|---|
| PALM22E | T5 | O_T | 2,656 | 2,671 | 1,214 | 2,671 |
| | | FDR | 0.11 ✓ | 0.017 ✓ | 0.0002 ✓ | 0.017 ✓ |
| | | efficiency | 0.80 | 0.74 | 0.34 | 0.74 |
| PALM22E | T4 | O_T | 4,386 | 2,530 | 1,217 | 2,514 |
| | | FDR | **0.98 ✗** | 0.065 ✓ | 0.004 ✓ | 0.057 ✓ |
| | | efficiency | 0.98 | 0.61 | 0.29 | 0.60 |
| PALM22E | T3 | O_T | 1,016 | 949 | 185 | 948 |
| | | FDR | 0.15 ✓ | 0.044 ✓ | 0.002 ✓ | 0.043 ✓ |
| | | efficiency | 0.79 | 0.63 | 0.12 | 0.63 |
| PALM09E | T1 | O_T | 1,272 | 34 | 0 | 34 |
| | | FDR | **0.89 ✗** | **0.25 ✗** | **∞ ✗** | **0.25 ✗** |
| | | efficiency | 0.88 | n.e. | n.e. | n.e. |

Mixture fits over covered tumour cells:

- PALM22E, v3 counts: T5 π̂ 1.00, f̂ 0.49; T4 π̂ 0.98, f̂ 0.18; T3 π̂ 0.95, f̂ 0.54.
- PALM09E T1, v3 counts: f̂ = 0.004 < ρ₀ + 0.05, so efficiency is **not evaluable**. The fit finds no MUT component.
- `main`'s high efficiencies at PALM22E T4 and PALM09E T1 come with ρ₀ = 0.49 and 0.20. The fit there separates two high-MUT populations, not mutant from wild type.

**3a**, v3 − (v2 + min2), needs ≥ 40 calls:

- PALM22E: T5 +1,457, T4 +1,297, T3 +763, all **pass**.
- PALM09E: T1 +34, a **tie** (3a not met).

**3c undetermined fraction**, as covered cells with an NA call (Wilson CIs in the TSVs):

| sample | target | compartment | `main` | v2 / v3 | v2 + min2 |
|---|---|---|---|---|---|
| PALM22E | T5 | Leukemic | 0.13 | 0.16 | 0.56 |
| PALM22E | T5 | Negatives | 0.50 | 0.63 | 0.64 |
| PALM22E | T4 | Leukemic | 0.01 | 0.10 | 0.41 |
| PALM22E | T3 | Leukemic | 0.30 | 0.34 | 0.83 |
| PALM09E | T1 | Leukemic | 0.02 | 0.12 | 0.14 |
| PALM09E | T3 | Leukemic | 0.57 | 0.81 | 0.86 |

## One-UMI test (QC cells; MUT fraction among cells with exactly one labelled UMI)

| sample | target | v3 counts: Neg | v3 counts: other | `main` counts: Neg | `main` counts: other |
|---|---|---|---|---|---|
| PALM22E | T3 | 0.008 (254) | 0.51 (1,121) | 0.051 (235) | 0.49 (1,017) |
| PALM22E | T4 | 0.003 (303) | 0.16 (541) | **0.58** (165) | 0.65 (169) |
| PALM22E | T5 | 0.005 (765) | 0.49 (1,166) | 0.016 (624) | 0.47 (939) |
| PALM09E | T1 | 0 (187) | **0** (196) | 0.085 (71) | 0.23 (62) |
| PALM09E | T2 | 0 (243) | 0.007 (289) | 0.005 (216) | 0 (337) |
| PALM09E | T3 | 0.080 (201) | 0.049 (226) | 0.16 (324) | 0.18 (430) |

Numbers in brackets are cell counts; Wilson CIs are in `one_umi_<S>.tsv`.

- At PALM22E T3–T5 a single MUT UMI in the v3 counts is informative: negatives read MUT ≤ 0.8 % against 16–51 % elsewhere.
- In `main`'s counts it is not informative at T4.
- At PALM09E T1 the v3 counts carry no MUT single UMIs at all.
- At PALM09E T3, negatives read MUT more often than the other cells, which fits the map's germline/artefact suspicion for this target.

## Deviations from the pre-registration

1. **`main` p check restricted to m ≥ 1.** Where m = 0, `main`'s `nbinom.cdf(k; 0, e)` is NaN, while the imported `call` gives p = 1 (k > 0). Neither is MUT, and both rules call WT iff k ≥ 2. The check therefore compares p on m ≥ 1 (max diff 1e-16) and separately requires `main`'s p to be NaN wherever m = 0. As an addition, it also requires the imported `call` with θᵢ to reproduce `main`'s `genotype` column exactly. Both passed.
2. **`main`'s mixture fit uses per-cell ρᵢ.** The pre-registration states the EM with a single ρ₀. For `main`, ρᵢ = θᵢ + c_s^main varies by cell. The likelihood uses ρᵢ per cell. The constraint f > ρ₀ and the "not evaluable" test f̂ < ρ₀ + 0.05 use the median ρᵢ over covered tumour cells (reported as `rho0`).
3. **`main`'s covered cells.** "Covered" uses each version's own counts, as pre-registered. `main`'s calibration set is the same B barcodes, restricted to those covered in `main`'s counts, so its B cell counts differ slightly from v3's. Evaluability uses v3's counts, as pre-registered.
4. **One-UMI test from both count sources.** The pre-registration names v3's diagnostic. It is given for the v3 run's counts and, for comparison, `main`'s. "Other" = all non-negative QC cells.
5. **Synthetic smoke data.** The sweep's (#68) synthetic generator re-seeded its random generator for every base, so its 3,000 barcodes collapsed to 4 distinct cells. Its smoke test was weaker than reported. Its real-data checks are unaffected. The generator is fixed here (`negative-control/mk.py`).

Items 6–8 were found by an audit during the verdict (#69). None changes a verdict.

6. **Clonal set fixed after the exploratory sweep.** The pre-registration requires the clonal set to be fixed before scoring starts, and it was. But #63 fixed it after the sweep's (#68) provisional verdicts were known, and it used v3 counts to move PALM22E T4 in and T7 out. That is a forking path. No verdict depends on it: `main` already fails PALM22E on criterion 1, and T7 would pass 2 and 3a for v3 (FDR 0.08, +176 calls over v2 + min2).
7. **Outlier flag at c_s = q_T = 0.** At PALM22E T2, c_s = q_T = 0, so the rule c_s ≥ 0.5 q_T is met literally, yet the table shows no flag. No effect: the target is neither evaluable nor clonal, and it has no tumour MUT calls.
8. **Calibration passes when no calls are made.** The one-sided calibration test passes any rule that makes no MUT calls in B (O_B = 0). v2 + min2 passes PALM22E calibration this way (at T4 in only 11 of 20 seeds). This is consistent with the pre-registered rule, but the earlier note that the ratio is "1 by convention" was wrong; it is 0.

## Reading for the verdict ticket

- Under the pre-registered bar, **v3 and v2 + min2 are the only contenders that pass anywhere (PALM22E)**, and **no contender passes PALM09E**.
- On PALM22E, v3 strictly dominates v2 + min2 in sensitivity: ×2.1 tumour calls at T5 and T4, ×5 at T3. It is also calibrated where v2 + min2 is calibrated only by abstaining.
- `main`'s failure is structural, not marginal:
  - its per-cell error is floored at 1e-4 in most cells;
  - its counts at the insertion targets read MUT in a large fraction of negative-cell UMIs.
- PALM09E's failure has two separable parts:
  - **T1, counts**: no rule over the v3 counts can pass. #70 found weak tumour signal rather than a v3 defect; T1 stays clonal and its fail stands (see [Verdict](#verdict)).
  - **T3, call rule**: v3 fails both calibration conditions, with 5 of 20 seeds passing. The exploratory v3-ub passes it in-sample.

## Verdict

Decided in *Verdict: which genotyping version to use for blast calls* (#69). The recommendations were reviewed by three independent reviewers (blind, adversarial, pre-registration audit). The audit re-derived all eight pass/fail outcomes from the pre-registered rules and found the table above faithful.

**Against the pre-registered bar, no contender passes both samples, so the map's destination is only partly met.** v3 passes PALM22E and fails PALM09E; v2 + min2 does the same. `main` and v2 fail both.

**Forward version (provisional): genotype-v3 with `--negative_cells` = T/B/NK, `min_mut_umis = 1`.**

- **Why v3.** It is the best contender: on PALM22E it is calibrated (held-out O/E 0.91 at T4, 0.97 at T5) and keeps 2–5× the tumour MUT calls of v2 + min2.
- **Why provisional.** It passed one sample out of two. Its calibration is shown only at low contamination (c_s 0.5–0.6 %, at PALM22E T4 and T5, passing in 17 and 14 of 20 seeds). The one higher-contamination evaluable target (PALM09E T3, c_s 4 %) fails.
- **Rejected.**
  - `main`: miscalibrated in negatives by 250× to 1,200× (held-out O/E 249–1,210). It must not be used at indel targets.
  - v2 (v3 without negative cells): miscalibrated by 6× to 400×. v3 must always be run with negative cells, and the sample name must match the negative-cell list, or c_s silently becomes 0.
  - v2 + min2: it calibrates only by making no calls in B (deviation 8), keeps 2–5× fewer tumour calls (efficiency 0.12–0.34), and still fails PALM09E T3 (O/E 2.4e4).

**Validated use: clone-fraction estimates from high-prevalence clonal calls, per target.** α = 0.1 per cell brings vagrant-RNA MUT calls down to the nominal rate; it does not remove them. In-sample false-MUT rates in negatives are 0.3–7.4 %. Tumour FDR is ≤ 0.06 at PALM22E's clonal targets. **Not validated:**

- multi-target "any MUT" cell labels. With three targets at 3–7 % each, roughly 10–20 % of WT cells would be labelled MUT somewhere;
- low-burden, MRD-like or T/B-involvement questions, where FDR grows as prevalence falls.

**PALM09E outcomes as recorded:**

- **T1, criteria 2 and 3a: fail**, as pre-registered. §5.2 makes O_T-driven failure a fail, and #63 kept T1 clonal knowing its signal was weak. *Interpretation, not a change of outcome:* #70 shows the fail reflects too little tumour signal in the counts (+4 points of tumour excess in reads, against +14 at PALM22E T4), not a v3 defect.
- **T3, criterion 1: fail** on both conditions (median p 0.048, median O/E 1.70; 5 of 20 seeds pass; in-sample false-MUT rate in negatives 7.4 %). **The cause is unresolved.** Two hypotheses stand:
  - B-cell over-dispersion: B-cell c is 13.5 % against T-cell 3.4 %, but that rests on 5 of 37 UMIs;
  - a germline or artefact target that the outlier flag did not catch: in the one-UMI test, negatives read MUT more often than other cells (0.080 against 0.049).
- **T2:** flagged and excluded. v3 passes it anyway (O/E 1.03), so no germline/clonal-haematopoiesis check is needed for this verdict.

**Follow-up, as one separate effort with a fresh pre-registration and new samples:**

- v3-ub, as a *candidate*: chosen on PALM09E T3 (O/E 1.70 → 1.29 there), it over-corrects elsewhere (O/E 0.27–0.69) and costs ≤ 3.5 % of tumour calls. It cannot be confirmed on the data it was chosen from.
- a germline/artefact check of PALM09E T3;
- position-tolerant indel matching, as its own ticket: a code change, measured only at read level, which adds noise at the 1 bp insertion target (PALM09E T2 negatives 27 → 93 cells).

## Reproduce

- **Scripts** (in the pinned v3 image, with `GENOTYPING_DIR` = `bin/long_reads` of `genotype-v3`@`61d8809`):
  - `negative-control/score.py`: the pre-registered scoring;
  - `negative-control/score_lib.py`: shared computations;
  - `negative-control/mk.py` and `mk_main.py`: synthetic smoke inputs; `main`'s output is made by `main`'s own `genotyping.py`.
- **Outputs** are in `negative-control/score-out/`:
  - `score_<S>.tsv`: every §6 column, per target × contender;
  - `one_umi_<S>.tsv`, `checks_<S>.json`, `verdict_<S>.json`.
- **HPC only:** input paths, the T-label key and the sbatch wrapper are in the working folder's `q67_score/` (`scripts/run.sbatch`).
