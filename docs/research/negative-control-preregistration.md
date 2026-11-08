# Pre-registration: T/B/NK negative controls vs blast calls (PALM22E, PALM09E)

Wayfinder ticket: *Pre-register the pass bar and split seeds* (#64), map *T/B negative controls vs blast calls* (#60).

This file is committed **before any full run or scoring**. It fixes the pass bar, the held-out split, the evaluability thresholds, the outlier rule and the scoring computations, so the verdict cannot be tuned to the results. Scoring (#67) cites this file's commit SHA. Any deviation must be listed in the report with its reason. The exploratory sweep (#68) is outside the verdict and is not bound by the decision rules below, but it uses the same computations.

**Patient data.** The repo is public. Everything here and in the report is aggregate. Targets are **T1…Tn** per sample, following the key in the HPC working folder's README. No barcodes, gene names, coordinates or patient-data paths appear in the repo.

## 1. Inputs and cell sets

- **Runs** (#66), from the canonical tagged BAMs (#62), with samplesheet sample names `PALM22_E` / `PALM09_E`:
  - one genotype-v3 run per sample, at `genotype-v3`@`61d8809`, with `--negative_cells` = the T/B/NK CSV and the default `min_mut_umis = 1`;
  - one `main` run per sample, at `main`@`9fe147a` or a later `main` commit recorded in #66.
- **Cell sets** per sample, all restricted to **QC cells** and matched on the bare 16-mer barcode:
  - **Negatives** (`Neg`) = QC ∩ T/B/NK list (= T/NK ∪ B of the annotation, #61).
  - **Tumour** (`Leuk`) = the annotation's Leukemic group.
  - **Ambiguous**: reported only, never scored against the bar.
  - Non-QC barcodes are dropped before anything is computed.
- **Covered.** A cell is covered at a target when its row has n = WT + MUT ≥ 1 *labelled* UMIs, counted from that version's own counts. OTH UMIs do not count towards n, because the call uses only m = MUT and k = WT. A cell with no row at a target is uncovered.
- **Clonal targets.** These come from *Clonal targets in PALM22E and PALM09E* (#63) and are fixed before scoring starts. Criteria 2 and 3 apply only to them.

## 2. Contenders and their call rules

Every contender uses the same one-sided binomial tail, with p = P(X ≥ m), X ~ Binomial(n, θ), and calls **MUT iff m ≥ m_min and p < 0.1**. This holds for `main` too: its `nbinom.cdf(k; m, e)` is identical to P(X ≥ m | Binomial(m + k, e)). The contenders differ only in the nominal null θ, in m_min and in which counts they use:

| contender | counts | nominal null θ | m_min |
|---|---|---|---|
| `main` (v1) | `main` run | θᵢ = max(errᵢ, baseline_error): per cell, the MIS-based `err` column floored as in the run | 1 (m = 0 is never MUT) |
| v2 | v3 run's consensus counts | P_s (target-pooled, OTH-based; from the v3 run's `_targets.csv`) | 1 |
| v2 + min2 | v3 run's consensus counts | P_s | 2 |
| v3 | v3 run's consensus counts | P_null = min(1, P_s + c), where c is c_A in the held-out split and c_s in-sample | 1 |

For a rule (θ, m_min) and n UMIs, the MUT region is m ≥ m*(n), where m*(n) is the larger of m_min and the smallest m with P(X ≥ m | Bin(n, θ)) < 0.1 (m* = ∞ if there is none). The rule's probability of calling a cell MUT when the true per-UMI null MUT rate is ρ is:

> **P(MUT | n; θ, m_min, ρ) = P(X ≥ m*(n)), X ~ Binomial(n, ρ)**  (0 if m*(n) > n)

- **Nominal false-MUT probability:** ρ = θ. This is the rule's own claim about itself, ≤ 0.1 per cell by construction.
- **Estimated true false-MUT probability:** ρ = P_s + c_s for v2, v2 + min2 and v3, and ρᵢ = θᵢ + c_s^main for `main`. Both are capped at 1. Here c_s is the in-sample contamination rate from all negatives, computed exactly as genotype-v3's `target_contamination`: Σ MUT / Σ (WT + MUT) over the negatives, with 0 if the denominator is 0. c_s^main is the same quantity from `main`'s counts.

**Implementation checks (abort on failure):**

- Re-calling the v3 run's counts with the in-sample P_null must reproduce its `genotype` column exactly.
- Recomputing `main`'s p from (m, k, θᵢ) must match its `p` column to 1e-9.
- The call function must be imported from `bin/long_reads/genotyping.py` at `genotype-v3`@`61d8809`, never re-implemented.

## 3. Held-out split

- **Split.** Per sample, the negatives are sorted lexicographically by barcode. For each seed: `rng = numpy.random.default_rng(seed)`, `perm = rng.permutation(len(Neg))`, **A** = `Neg[perm[: len(Neg) // 2]]`, **B** = the rest. The split is unstratified and at the cell level, and every target shares it. The scoring container's numpy version is recorded in the report.
- **c_A** = Σ_A MUT / Σ_A (WT + MUT) at the target, from the v3 run's counts, with 0 if the denominator is 0. v3 cells in B are re-called with P_null = min(1, P_s + c_A). P_s comes from the run and is not re-estimated: it is OTH-based, so it does not leak MUT information from B. v2, v2 + min2 and `main` do not use c, so their calls in B are the same for every seed. They are still scored on the same B cells.
- **The 20 seeds** (drawn once with `random.Random(20261108).sample(range(1, 2**31), 20)`):

```
2119795751  703054139  967594051 2137092682  479082201
1922421208  640385892  828295203 1190487198   31658864
1178418751 1637046186 1438843744 2114041616  913149788
 174283017 1834471605  912037160 1502725731 2038317252
```

- **Across seeds** the report gives the **median and 5–95 % range** of every held-out quantity. With 20 seeds that is the mean of the 10th and 11th order statistics, and the 1st and 19th order statistics (numpy `percentile` with the default linear method). Pass/fail is decided on the medians (§5). The fraction of seeds that pass is reported next to it.
- **In-sample** (c_s from all negatives, scored on all negatives, i.e. what the v3 run itself produced) is reported alongside, but it does not decide the bar.

## 4. Evaluability and outliers

- **Evaluable seed:** at the target, ≥ 100 covered cells in B **and** ≥ 200 WT + MUT UMIs in A, counted from the v3 run's counts.
- **Evaluable target:** evaluable in ≥ 10 of the 20 seeds. Held-out medians are taken over the evaluable seeds only. A target that is not evaluable is reported as **not evaluable** (not failed) and drops out of criterion 1.
- **Outlier c_s (flag):** a target is flagged when its in-sample c_s from the v3 run is ≥ 0.25, **or** ≥ 0.5 × q_T. Here q_T is the MUT fraction of WT + MUT UMIs pooled over covered tumour cells, so the second test catches a target that does not discriminate negatives from tumour. Flagged targets are still scored and reported, but they are **excluded from the verdict** (germline or clonal haematopoiesis is not checked here; see the map's fog).

## 5. Pass bar (per contender, per target)

Write O for the observed MUT calls, and E for the sum of P(MUT | nᵢ; …) over the stated cells.

1. **Calibration (specificity)**, every evaluable, unflagged target. Computed per seed on the covered cells of B:
   - O_B = MUT calls in B.
   - E_B = Σ_B P(MUT | nᵢ; θ, m_min, ρ = θ), using the nominal null with c_A for v3.
   - p_cal = P(S ≥ O_B), where S is the Poisson-binomial sum of independent Bernoulli(P(MUT | nᵢ)). This is the exact generalisation of the one-sided binomial test, computed by dynamic programming.
   - **Pass if median p_cal ≥ 0.05 and median O_B / E_B ≤ 1.5.** If E_B = 0 and O_B > 0, the ratio is ∞. If O_B = E_B = 0, the ratio is 1.
2. **Precision in the tumour compartment**, every clonal, unflagged target, over covered tumour cells and in-sample:
   - Conservative FDR = E_T / O_T, where E_T = Σ_Leuk P(MUT | nᵢ; θ, m_min, ρ = estimated true rate) (§2).
   - It is conservative because it treats every covered tumour cell as WT (π = 0).
   - **Pass if FDR ≤ 0.2.** If O_T = 0, FDR = ∞ (fail).
3. **Sensitivity (secondary)**, clonal, unflagged targets:
   - **3a (gating, v3 only):** O_T(v3) − O_T(v2 + min2) ≥ 40 calls, the PALM14E run-to-run noise floor applied per target. A difference below 40 is a **tie**: 3a is not met and the report marks it as a tie. 3a does not apply to the other contenders.
   - **3b (reported against the bar, not gating):** detection efficiency ≥ 0.7, where efficiency = (O_T − (1 − π̂) · E_T) / (π̂ · n_cov,T).
     - π̂ and f̂ are the maximum-likelihood fit, over covered tumour cells, of mᵢ ~ (1 − π) Bin(nᵢ, ρ₀) + π Bin(nᵢ, f). Here ρ₀ is the estimated true null rate, held fixed, with f ∈ (ρ₀, 1] and π ∈ [0, 1].
     - The fit is run by EM from π = 0.5, f = 0.5, until the change in log-likelihood is < 1e-8, with at most 10,000 iterations.
     - It uses the v3 run's counts for v2, v2 + min2 and v3, and `main`'s counts for `main`.
     - Efficiency is not evaluable if π̂ < 0.05 or f̂ < ρ₀ + 0.05. Undetermined calls count as not detected.
   - **3c (reported):** the undetermined fraction (covered cells with an NA call ÷ covered cells), for every contender, compartment and target.

**Sample verdict per contender:** the contender **passes** on a sample when:

- every evaluable, unflagged target passes 1;
- every clonal, unflagged target passes 2;
- for v3, every clonal, unflagged target also passes 3a.

The ~14 calibration tests are not corrected for multiplicity: that leans towards failing, which is the safe side for a filter.

**Ties and replicates:** if a 3a tie, or a contender-vs-contender difference in O_T below 40 calls, decides which contender goes forward, the verdict ticket considers replicate runs (map fog). Otherwise there are no replicate runs.

**From pass bar to decision:** the pass bar only decides pass / fail / not evaluable per contender. Choosing among the contenders that pass is the verdict ticket's job (#69, HITL).

## 6. Reported quantities (all aggregate, per sample × target × contender)

- Coverage: covered cells in Neg, A, B, Leuk and Ambiguous; WT + MUT UMIs in A.
- c_s (in-sample), the seed-median c_A and its 5–95 % range, P_s, θ summary for `main` (median θᵢ), q_T, and the outlier flag.
- Calibration: O_B, E_B, O_B / E_B, p_cal (median, 5–95 %); the observed false-MUT rate O_B / n_cov,B with a **Wilson 95 % CI**; the same in-sample over all negatives.
- Precision: O_T, E_T and the conservative FDR. E_T with ρ = θ (nominal) is also given for comparison.
- Sensitivity: O_T for every contender, the v3 − (v2 + min2) difference, π̂, f̂, detection efficiency, and undetermined fractions with Wilson 95 % CIs.
- One-UMI test (`one_umi_mut_frac_neg`, `one_umi_mut_frac_other`), recomputed on QC cells, with Wilson 95 % CIs.
- Per-target outcome per criterion: pass / fail / tie / not evaluable / flagged. Then the sample verdict per contender.

Wilson intervals use z = 1.96. Ratios with a zero denominator are reported as ∞ or undefined, never dropped.
