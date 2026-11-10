# Pre-registration: PALM09E T3 cause (germline, artefact or B-cell over-dispersion)

Wayfinder ticket: *PALM09E T3 cause: germline, artefact or B-cell over-dispersion* (#73), map *Confirm genotype-v3 vs v3-ub on new diagnosis samples* (#72).

This file is committed **before any new look at the data**. It fixes the measurements, the thresholds, the decision rule and what each outcome means for the confirmation. The report cites this file's commit SHA and lists any deviation with its reason.

**Already known when this was written** (published aggregates from the scoring and sweep reports, #67 and #68): at T3, c_s = 4.0 % (32 MUT of 808 negative UMIs), and the MUT rate by negative group is T 3.4 % (24/698), NK 4.1 % (3/73) and B 13.5 % (5/37, 20 cells). q_T = 13 %, and the tumour mixture fit puts π̂ ≈ 2.9 % of covered tumour cells at f̂ ≈ 0.98. In-sample, 23 negatives are called MUT against 13.9 expected; held out, the median O/E is 1.70 and the median p is 0.048. The one-UMI test reads MUT in 8.0 % of negatives and 4.9 % of other cells. T3 is a genomic G>A SNV. Nothing else about T3 has been looked at.

**Patient data.** Aggregate numbers only: targets as T1…Tn, no barcodes, gene names, coordinates or patient-data paths.

## 1. Inputs and cell sets

- **Counts:** the genotype-v3 full run of PALM09E (#66, `genotype-v3`@`61d8809`, T/B/NK negatives, `min_mut_umis = 1`), `<S>.csv`.
- **Reads:** the canonical tagged BAM of PALM09E (#62). Primary alignments with MAPQ ≥ 20, carrying a CB tag, whose base at the T3 position is REF or ALT. Deletions and other bases at the position are counted but take no part in the tests.
- **Cell sets.** All restricted to QC cells and matched on the bare 16-mer barcode, as in the scoring (#67):
  - **T** = negatives with annotation `aggregated_ct` starting with `T_`;
  - **NK** = negatives with `aggregated_ct` = `NK_cells`;
  - **B** = negatives with `ct` = `B_cells`;
  - **Neg** = T ∪ NK ∪ B (= the T/B/NK list ∩ QC);
  - **Leuk** = `ct` = `Leukemic`; **Amb** = `ct` = `Ambiguous`; **Other** = the remaining QC cells (reported only).
- **Covered:** a cell is covered at a target when n = WT + MUT ≥ 1 in the v3 counts.

## 2. Measurements

**Count level (v3 counts, T3):**

- **R (lineage rates).** For each of T, NK, B, T/NK, Neg, Leuk, Amb and Other: the covered cells, UMIs, MUT UMIs, the UMI MUT rate r = Σ MUT / Σ (WT + MUT) with a Wilson 95 % CI, and the fraction of covered cells with ≥ 1 MUT UMI.
- **D1 (over-dispersion in negatives).** Fit, over covered Neg cells, a binomial model (mᵢ ~ Bin(nᵢ, r)) and a beta-binomial model (mᵢ ~ BetaBin(nᵢ, r, ρ), ρ ∈ [0, 1)) by maximum likelihood. LR = 2 (ℓ_BB − ℓ_Bin), and p = ½ P(χ²₁ ≥ LR), the boundary-corrected test. Report ρ̂. The same fit is reported for T alone and for Leuk (descriptive).
- **D2 (B-specificity).** One-sided Fisher exact test (B greater) on covered cells with ≥ 1 MUT UMI vs without, B vs T/NK.
- **U (lineage homogeneity).** Pearson χ² test of homogeneity of the UMI MUT rate across T/NK, B and Leuk (2 × 3 table of MUT vs WT UMIs), with a Monte Carlo p-value (10⁵ draws, seed 73) when any expected count is < 5.

**Read level (tagged BAM, T3 position):**

- **A1 (strand).** In Neg reads: the 2 × 2 table of allele (ALT, REF) × aligned strand (+, −). Fisher exact two-sided p, and Δ = f₊(ALT) − f₊(REF).
- **A2 (read position).** In Neg reads: the distance from the T3 base to the nearer end of the aligned part of the read, in query bases. ALT vs REF: Mann–Whitney two-sided p, and the ratio of medians ALT / REF.
- **A3 (molecule purity).** Raw molecules are (CB, UB) pairs with ≥ 3 qualifying reads. Among molecules whose ALT read fraction is ≥ 0.5, compare the ALT fraction in Neg with Leuk: one-sided Mann–Whitney (Neg lower) and the difference of medians, Leuk − Neg.
- **A4 (context).** In Neg reads: the per-read rate of the same substitution (reference G read as A) at T3 vs at every other reference-G position within ±100 bp of T3 that has ≥ 200 qualifying Neg reads. Report the 95th percentile of the background and T3's ratio to it, separately for all G positions and for G positions sharing T3's trinucleotide context.

## 3. Signals

| signal | holds when |
|---|---|
| **G** (germline) | for each of T, B and Leuk (and NK if it has ≥ 20 covered cells), the Wilson CI of r overlaps [0.3, 0.7] (heterozygous); or every lower bound is ≥ 0.9 (homozygous) |
| **A** (artefact) | any of: **A1** p < 0.01 and \|Δ\| ≥ 0.25; **A2** p < 0.01 and ALT/REF median ≤ 0.5; **A3** p < 0.01 and difference ≥ 0.15; **A4** T3's rate ≤ 3 × the 95th percentile of the all-G background |
| **D** (over-dispersion) | **D1** p < 0.05, or **D2** p < 0.05 |
| **U** (uniform) | U's p ≥ 0.05 |

## 4. Decision rule

Applied in order; the first row that matches decides.

| # | condition | cause |
|---|---|---|
| 1 | G | **germline** |
| 2 | A and not D | **artefact** |
| 3 | U, not A and not D | **uniform signal**: germline with allelic imbalance, or a molecule-level site artefact (not separable on RNA). Same consequence as germline or artefact |
| 4 | D and not A | **over-dispersion**; labelled **B-specific** if D2 holds, **cell-level, not B-specific** if only D1 holds |
| 5 | anything else (including A and D together, or no signal) | **inconclusive** |

**Power, stated in advance.** B has 20 covered cells at T3, so D2 detects only a large B excess. A D2 that fails is not evidence against B-specificity.

## 5. Consequences for the confirmation (map #72)

- **Over-dispersion:** v3-ub, or a per-lineage or beta-binomial null, stays a contender in the confirmation pass bar.
- **Germline, artefact or uniform signal:** the outlier flag is under-sensitive at T3, and v3-ub was tuned to a non-problem. The fog patch *Outlier-flag revision* graduates; v3-ub may still be scored on the new samples, but T3 gives it no support.
- **Inconclusive:** v3 and v3-ub both go to the new samples unchanged, and the PALM09R pair check graduates from the fog.

## 6. Reported, not deciding

- Cross-target: for every other PALM09E target with ≥ 20 covered B cells, r(B) / r(T/NK). Vagrant contamination that favours B cells predicts a ratio > 1 there too.
- MUT-carrying vs other covered negatives: total UMIs over all targets, and the fraction flagged as doublets by the annotation (`scDblFinder.class`).
- Base quality of ALT vs REF bases (medians); the A1–A3 measurements in Leuk reads, for reference.
- Editing: on the transcript, a genomic G>A on a minus-strand gene reads as C>U, an APOBEC-type change. RNA data cannot separate C>U editing from a somatic or germline change. The report states the transcript strand and, where it applies, names editing as a caveat on any non-artefact outcome.

## 7. Implementation

- One Python script, `docs/research/t3-cause/t3_cause.py`, run on `hpc_new` by sbatch in the pinned NanoGen image used by the scoring (#67). Output is aggregate only and goes to `docs/research/t3-cause/out/`.
- The script is checked first on a synthetic fixture with a known over-dispersed negative group, which must give a D1 p < 0.05.
- Implementation check (abort on failure): the Neg pooled rate recomputed from the counts must equal the run's c_s at T3 to 1e-12.
