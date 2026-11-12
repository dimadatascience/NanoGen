# Pre-registration: confirming genotype-v3 vs v3-ub on new diagnosis samples

Wayfinder ticket: [Pre-register the confirmation pass bar](https://github.com/dimadatascience/NanoGen/issues/78), map [Wayfinder: confirm genotype-v3 vs v3-ub on new diagnosis samples](https://github.com/dimadatascience/NanoGen/issues/72). Compare the earlier pre-registration, [`negative-control-preregistration.md`](https://github.com/dimadatascience/NanoGen/blob/76c3f99/docs/research/negative-control-preregistration.md) (`76c3f99`), cited below as **P64**.

This file is committed **before any genotype run or scoring on the census**. It fixes the contenders, the cell and target sets, the held-out split and its seeds, the pass bar, the decision rule and every reported quantity, so the verdict cannot be tuned to the results. Scoring cites this file's commit SHA and lists every deviation with its reason.

**Patient data.** The repo is public. Everything here and in the report is aggregate. Targets are **T1…Tn** per patient, following the key in the HPC working folder's README. No barcodes, gene names, coordinates or patient-data paths appear in the repo.

## 1. Claim and its stated limits

The verdict decides whether **genotype-v3** holds as NanoGen's forward setting for **per-target clone-fraction estimates from high-prevalence clonal calls**, and whether **v3-ub** replaces it.

- The calibration half of the verdict is a **pooled calibration test**. A pass means "no gross *average* miscalibration on new patients". **Per-target calibration is not confirmable on this census and stays provisional.**
- **Blind spot (reported verbatim in the verdict).** On this census, neither the pooled test nor the site screen can rule out a single-site miss the size of PALM09E T3's (O/E ≈ 1.7). A pooled pass therefore means "no average miss ≥ 1.5× and no gross site miss on well-covered targets", **not** "no T3-like site misses". See [`pool-guard.md`](https://github.com/dimadatascience/NanoGen/blob/cc91417/docs/research/pool-guard.md): no pool target reaches E_B 5.7 at c = 1 %, the screen's 80 % power point for a 3× miss.
- The tumour-side half covers only the tumour-evaluable loci of §6.1, and the verdict says which.

## 2. Inputs and cell sets

- **Patients (frozen pool, 8):** PALM07E, PALM08E, PALM11E, PALM14E, PALM21E (census) and AML4, AML5, sAML1 (old-panel FASTQ restart).
- **BAMs:** the final canonical tagged BAMs. These are the regenerated short-read-known-list BAMs for PALM08E, 11E, 14E and 21E ([#85](https://github.com/dimadatascience/NanoGen/issues/85)), library 2 for PALM07E ([#88](https://github.com/dimadatascience/NanoGen/issues/88)), and the FASTQ-restart BAMs for AML4, AML5 and sAML1 ([#84](https://github.com/dimadatascience/NanoGen/issues/84)). Paths, QC lists, negative lists and tumour labels are in the HPC README sections for those tickets.
- **Cell sets,** per patient, all restricted to **QC cells** and matched on the bare 16-mer barcode:
  - **Negatives (`Neg`)** are the T/B/NK cells from each intake.
  - **Tumour (`Leuk`)** is the intake's Leukemic label. For AML4, AML5 and sAML1 this is the older study's phenotype-only tier-I call.
  - **Ambiguous** cells are reported only.
- **Covered** is as in P64 §1: n = WT + MUT ≥ 1 labelled UMIs, from that run's own counts.
- **Late inputs are frozen out.** No clonal set for PALM07E or PALM21E, and PALM05E stays out of the census. If bulk VAF or an annotation arrives later, it is scored only in the supplementary arm (§9).

## 3. Runs and contenders

**Runs,** one per patient each, on `hpc_new` via sbatch + singularity. Mock runs size resources, and full runs request ≤ measured + 50 %.

- **v3 run:** `genotype-v3`@`61d8809`, `--negative_cells` = that patient's negatives, `min_mut_umis = 1`.
- **Indel-arm run:** `prototype/position-tolerant-indel`@`c77e147`, `--indel_window 12`, otherwise identical. This branch is never merged.

**Contenders.** All use the P64 §2 call rule: MUT iff m ≥ m_min and p < 0.1, with p = P(X ≥ m), X ~ Bin(n, θ), and m_min = 1. Here c = Σ MUT / Σ (WT + MUT) over a set of negatives, 0 if the denominator is 0, and CPU(x, n) = the one-sided 95 % Clopper–Pearson upper bound = `scipy.stats.beta.ppf(0.95, x + 1, n − x)`, with 1 if x = n or n = 0.

| contender | counts | nominal null θ (held-out) | nominal null θ (in-sample) |
|---|---|---|---|
| **v3** (incumbent) | v3 run | min(1, P_s + c_A) | min(1, P_s + c_s) |
| **v3-ub** (challenger) | v3 run (a re-call, so it has no run noise relative to v3) | min(1, P_s + CPU(x_A, n_A)) | min(1, P_s + CPU(x_s, n_s)) |
| **indel prototype** (secondary arm, §8) | indel-arm run | the rule chosen by §7, applied to its own counts | same |

**Implementation checks (abort on failure):**
- the call function is imported from `bin/long_reads/genotyping.py` at `genotype-v3`@`61d8809`, never re-implemented;
- re-calling each run's counts with its in-sample P_null reproduces that run's `genotype` column exactly;
- the scoring library is frozen at a SHA recorded in the report.

**Estimated true null** (for FDR, both contenders): ρ = min(1, P_s + c_s) with point-estimate c_s, as in P64 §2.

## 4. Held-out split and seeds

- **Split:** as in P64 §3. Negatives are sorted lexicographically by barcode, `perm = numpy.random.default_rng(seed).permutation(len(Neg))`, A = the first ⌊len/2⌋, B = the rest. The split is unstratified, at cell level, shared by all targets of a patient, and independent per patient.
- **20 new seeds**, drawn once with `random.Random(20261112).sample(range(1, 2**31), 20)`:

```
 452483223 1695653105 1085084613  441224510  821422349
1758734097 1576893678 1902659721  536475717 2055609293
1297200436 2099741074  248750689 1900119910 1632811357
 522860322  934946355  760461626  377745768 1704488278
```

- **Across seeds:** the median (mean of the 10th and 11th order statistics) and the 5–95 % range (numpy `percentile`, linear). Pass/fail is decided on the medians, and the fraction of passing seeds is reported beside each.

## 5. Calibration (primary, gating)

**Pool.** The 57 targets frozen by [Final pool](https://github.com/dimadatascience/NanoGen/issues/86) (`546386e`, `final-pool/out/pool86.tsv`, `in_pool = True`):

| patient | pool targets |
|---|---|
| PALM07E | T1–T5 |
| PALM08E | T1, T3–T7 |
| PALM11E | T1, T3–T5 |
| PALM14E | T1–T7 |
| PALM21E | T1–T5, T7–T12 |
| AML4 | T1–T4, T6, T9, T10 |
| AML5 | T1–T4, T7, T10, T11, T16 |
| sAML1 | T1, T2, T4, T7–T12 |

**Rules for building the pool:**
- Targets outside this list (proxy E_B = 0) are excluded even if they realise E > 0.
- Targets flagged by the outlier flag (§5.4) are excluded from the sums and reported.
- **Zero-information target-seeds:** a target-seed with n_A = Σ_A (WT + MUT) = 0 is left out of the sums **for both contenders**, both O and E. The count of such target-seeds is reported per patient. If it exceeds 10 % of the pool's target-seeds, the verdict says so.
- Otherwise, a member contributes its O_B even if its E_B is 0.

### 5.1 Pooled calibration test

For each seed and each contender, over every covered B cell at every included target-seed:
- O = Σ MUT calls;
- E = Σ P(MUT | nᵢ; θ, m_min, ρ = θ), using the contender's nominal held-out null;
- p = P(S ≥ O), where S is the exact Poisson-binomial sum, computed by dynamic programming.

**Pass if the seed-median p ≥ 0.05 and the seed-median O/E ≤ 1.5.** Ratio conventions follow P64 §5: E = 0 with O > 0 gives ∞, and O = E = 0 gives 1.

- **No lower bound gates.** The lower-tail p = P(S ≤ O) and O/E are reported for both contenders. Over-correction is priced by the tumour-call loss cap (§6.4).
- **Non-gating fragility rerun:** the same test, also excluding target-seeds with n_A < 11 (the smallest n at which CPU(0, n) ≤ 0.25). If either contender's pass/fail differs between the primary and the rerun, the verdict is labelled **fragile**.

### 5.2 Per-patient reporting (non-gating)

For each patient and contender, report:
- the seed-median Σ O / Σ E, with an exact 90 % Poisson (Garwood) CI on O divided by E;
- the leave-one-patient-out pooled O/E and p.

### 5.3 Site screen (flags, never vetoes)

- **Per target, per seed:** the one-sided **exact** (not mid-p) Poisson tail p = P(Poisson(E_B) ≥ O_B). E_B = 0 with O_B > 0 gives p = 0.
- **Per target:** the median over that target's included seeds. A target included in no seed is reported as "no negative data" and is not tested.
- **Across targets:** Holm across the tested targets at α = 0.10, separately for each contender.
- **A flagged target is a site exception.**
  - It stays in Σ O / Σ E and does not fail the verdict.
  - The verdict lists it, with its O/E and both contenders' status.
  - It opens a follow-up germline/artefact check (the T3 pattern), recorded in the map's fog.
  - It does not change the decision rule (§7).
- **No post-hoc drops and no heterogeneity veto.**

### 5.4 Outlier flag (exclusion, pre-registered)

The flag is unchanged from P64 §4: in-sample c_s ≥ 0.25, **or** c_s ≥ 0.5 × q_T, both from the v3 run. Here q_T is the MUT fraction of WT + MUT UMIs pooled over covered tumour cells.
- **Edge case fixed:** c_s = q_T = 0 is **unflagged**. This was P64 deviation 7.
- Flagged targets are excluded from the pooled sums and from every tumour criterion, and are reported.
- **Molecule purity** (negatives' ALT in partial, not whole, molecules) is reported per target as descriptive only.

## 6. Tumour side (per target, gating)

### 6.1 Tumour-evaluable loci (fixed now, from the frozen proxy table)

A WES-clonal target ([#77](https://github.com/dimadatascience/NanoGen/issues/77), [#87](https://github.com/dimadatascience/NanoGen/issues/87): bulk WES CCF ≥ 0.8) is **tumour-evaluable** if it has **≥ 100 covered Leukemic cells** in `pool86.tsv` @ `546386e` (column `leuk_covered`). That is P64 §4's 100-cell minimum. A threshold of 50 would select the same set: no clonal target has 48–107 covered tumour cells.

**Tumour-evaluable: 17 targets = 15 loci in 5 patients.**

| patient | tumour-evaluable targets |
|---|---|
| PALM08E | T3, T4 |
| PALM11E | T1 |
| PALM14E | T1, T2, {T3–T5} (one phased locus), T6, T7 |
| AML5 | T1, T3 |
| sAML1 | T7, T9, T10, T11, T12 |

- **AML4 is calibration-only** (≤ 19 covered tumour cells at every clonal target), as are PALM07E and PALM21E (no WES).
- Every other clonal target is listed as **under-covered, not evaluable**. Under-covered targets never count toward the cap (§6.5).
- If a target passes here but has < 100 covered tumour cells in the v3 run, it stays in and the drop is reported.
- **Phased locus:** PALM14E T3–T5 is gated per target but **counted once** in every tally, and its O_T enters the pooled loss as the median of the three.
- **Thin null:** PALM14E {T3–T5} has 6–9 negative-control UMIs. It stays in the gate and the loss, and is marked "thin null" in the scope line (§6.3).

### 6.2 No-signal targets (frozen before any contender is scored)

- **Test:** on each tumour-evaluable, unflagged target, fit the P64 §5 3b mixture once on the **v3 run's counts**, with ρ₀ = P_s + c_s (point estimate): mᵢ ~ (1 − π) Bin(nᵢ, ρ₀) + π Bin(nᵢ, f). EM starts from π = f = 0.5 and stops at Δ log-likelihood < 1e-8, with at most 10,000 iterations.
- **No signal:** a target is no-signal if **π̂ < 0.05 or f̂ < ρ₀ + 0.05**. Both clauses are needed: on PALM09E T1 only the f̂ clause fires.
- The list is committed to the scoring branch before either contender's tumour outcomes are computed.
- No-signal targets are reported as **counts-level failures**, not caller failures, and leave the gate for both contenders.

### 6.3 Precision gate (both contenders)

On every tumour-evaluable, unflagged, signal-bearing target, in-sample over covered tumour cells:
- **Conservative FDR** = E_T / O_T, where E_T = Σ_Leuk P(MUT | nᵢ; θ, m_min, ρ = estimated true null).
- **Pass if FDR ≤ 0.2.** O_T = 0 gives ∞ (fail).
- **n_s = 0** (no negative-control UMIs in-sample): "null not estimable" for both contenders and out of the gate. This is moot on this census.
- **Weak-null check:** at the thin-null locus PALM14E {T3–T5}, v3's FDR is also reported with c_s replaced by CPU(x_s, n_s). If v3 passes only under the point estimate, the verdict carries a **weak-null pass** flag for that locus.

### 6.4 v3-ub tumour-call loss (replaces P64 criterion 3a)

- **Loss** = 1 − Σ O_T(v3-ub) / Σ O_T(v3), in-sample, over the **same gated set** as §6.3 (phased locus once). **Gate: loss ≤ 5 %.**
- The difference is exact: the same counts are re-called, so no tie band applies.
- Per-patient and per-target losses are reported. v3-ub's silence on targets outside the gated set is reported separately, as a deployment cost.
- v3 has no sensitivity gate of its own, because P64 3a's comparator (v2 + min2) is rejected.
- **Reported, not gating:** detection efficiency (P64 3b) and undetermined fractions (P64 3c).

### 6.5 Tumour side not evaluable

The tumour side is **not evaluable** if either:
- **3 or more of the 15 loci are no-signal**; or
- **fewer than 3 patients keep a gated locus.**

### 6.6 Scope line (mandatory in the verdict)

"Per-target precision confirmed on N covered, signal-bearing loci in K patients. U under-covered and M no-signal targets are listed. Thin-null loci: PALM14E {T3–T5}." Add a weak-null pass flag where it applies.

## 7. Decision rule

A contender **passes** when it passes the pooled calibration test (§5.1), every gated target passes §6.3, and the tumour side is evaluable (§6.5).

1. **v3 holds** (forward setting confirmed for the claim in §1) if v3 passes.
2. **v3-ub replaces v3** only if v3 **fails** the pooled calibration, **and** v3-ub passes, **and** v3-ub's loss ≤ 5 % (§6.4).
3. **If both pass, v3 holds.** The census cannot show the T3-scale benefit v3-ub was built for.
4. **Otherwise** (both fail, or v3 fails on the tumour side, or the tumour side is not evaluable), **neither is confirmed and v3 stays provisional.**

**Labels attached to the verdict, none of which changes it:**
- fragile (§5.1);
- site exceptions (§5.3);
- weak-null pass (§6.3);
- the blind-spot statement (§1);
- the scope line (§6.6).

## 8. Secondary arm: position-tolerant indel matching

The secondary arm is scored **only after the primary verdict is written**. It cannot rescue or sink the primary verdict.

- **Contender and baseline.** The contender is the indel-arm run (§3), called with the rule chosen by §7 (v3, or v3-ub if it replaced v3), with c_A, c_s and calls from its own counts. The baseline is the same rule on the v3 run.
- **"Worth a code ticket" only if all three hold:**
  1. **Gain.** Over the tumour-evaluable, signal-bearing clonal indel targets (PALM08E T3, PALM11E T1, PALM14E T3, sAML1 T11), the pooled in-sample tumour MUT calls rise by more than the summed per-sample tie bands **and** by ≥ 10 %.
  2. **No calibration harm.** The prototype passes §5.1 over the full pool. It has **no site exception at an indel target that v3 lacks**. The pool's indel targets are PALM08E T3, T6, T7; PALM11E T1; PALM14E T3; PALM21E T1, T2, T7; AML5 T2; sAML1 T11, of which the 1 bp insertions are PALM08E T6, PALM21E T2 and AML5 T2. The 10-target indel sub-pool O/E is reported, not gating.
  3. **Precision.** FDR ≤ 0.2 at every gated clonal indel target.

  Any failure gives "not worth a code ticket", with the reason.
- **Tie band (cross-run comparisons only):** band_s = max(⌈40 × rows_s / 34,979⌉, measured_s).
  - rows_s is the patient's cell × target row count in the v3 run. 40 calls in 34,979 rows is the PALM14E run-to-run floor.
  - measured_s is the number of SNV-target call differences between the two runs, scaled to all rows (× rows_s / SNV rows_s). At SNV targets the prototype is v3, so these differences are pure run noise.
  - A gain within the summed bands is a **tie**: not worth a code ticket, and **no replicate runs**.

## 9. Supplementary arm (late inputs)

If bulk VAF for PALM07E or PALM21E, or a Leukemic/Ambiguous annotation for PALM05E, arrives after this commit, it is scored under the identical rules of §5–§6 in a separate section. That includes the 100-cell coverage screen, read from the same proxy, re-run for PALM05E. These results are reported beside the verdict but **cannot change it**. A disagreement triggers a fresh confirmation effort.

## 10. Reported quantities (aggregate, per patient × target × contender)

- **Coverage:** covered cells in Neg, A, B, Leuk and Ambiguous; WT + MUT UMIs in A and in all negatives.
- **Contamination and the outlier flag:** c_s, CPU(x_s, n_s), the seed-median c_A and CPU(x_A, n_A) with 5–95 % ranges, P_s, q_T, and the outlier flag.
- **Calibration:**
  - O_B, E_B, O_B / E_B, upper- and lower-tail p (medians and 5–95 %);
  - the false-MUT rate with a Wilson 95 % CI;
  - site-screen p and Holm-adjusted outcome;
  - the excluded zero-information target-seeds.
- **Tumour:** O_T, E_T, FDR, π̂, f̂, efficiency, undetermined fractions with Wilson CIs, no-signal status, and v3-ub loss.
- **Indel arm:** tol_reads, tol_unseen_reads and brc_mismatch_umis summaries; calls per target in both runs; tie bands.
- **Outcomes:** pass / fail / tie / not evaluable / under-covered / no-signal / flagged / site exception per criterion, then the verdict and its labels.

Wilson intervals use z = 1.96.

## 11. Changes from P64, and why

| P64 | here | reason |
|---|---|---|
| per-target evaluability (≥ 100 B cells, ≥ 200 A UMIs) gates calibration | pooled test over 57 targets; site screen flags | no target is powered ([#76](https://github.com/dimadatascience/NanoGen/issues/76), [#83](https://github.com/dimadatascience/NanoGen/issues/83)) |
| flagged targets excluded; c_s = q_T = 0 ambiguous | same flag; edge case unflagged | P64 deviation 7 |
| 3a: v3 − (v2 + min2) ≥ 40 calls | v3-ub vs v3 loss ≤ 5 % | v2 + min2 rejected; the difference is exact |
| one-sided calibration only | loss cap prices over-correction; lower tail reported | P64 deviation 8 (calibration passes when no calls are made) |
| FDR gate on every clonal target | coverage screen, then no-signal, then per-target gate plus cap | PALM09E T1 (a counts failure) and 12 of 41 clonal targets with no covered tumour cells |
| absolute 40-call tie band | scaled per sample plus measured; indel arm only | primary contenders share counts |
| clonal set fixed after the sweep | clonal sets and tumour-evaluable list fixed here | P64 deviation 6 |
