# Clonal targets in the old-panel restart samples (AML4, AML5, sAML1)

Ticket: "Clonal targets in the old-panel restart samples" (issue 87), map "Wayfinder: confirm genotype-v3 vs v3-ub on new diagnosis samples" (issue 72).
Aggregate numbers only; targets are T1…Tn in panel-BED order (key in the HPC working-folder README).

## Source found

The [FASTQ restart power check](https://github.com/dimadatascience/NanoGen/issues/84) had concluded that no bulk WES existed for these patients. It does: an older DRAGEN tumour/normal WES batch (2023, hg38, same pipeline family as the PALM production batches) holds a diagnosis tumour (`_D`) and a matched normal (`_N`) for all three patients. The older single-cell study already used it for its WES-vs-single-cell clone-fraction comparison. The path is in the HPC README, not here.

The old panels were designed from this WES. Every panel target matches a WES call at the same position and allele, so the WES and the single-cell library are taken to be the same diagnosis sample. That is the same assumption the [new-diagnosis clonal ticket](https://github.com/dimadatascience/NanoGen/issues/77) made.

## Rule (unchanged from the new-diagnosis clonal ticket)

- PASS somatic call in the hard-filtered VCF, absent from the matched normal;
- purity p = 2 × median VAF of the top VAF cluster of copy-neutral autosomal panel targets. The top cluster is everything above the largest VAF gap;
- CCF = 2·VAF / p; **clonal ⇔ CCF ≥ 0.8** (point estimate);
- the clonal set is fixed before any contender is run on these samples. No v3, v3-ub or indel-prototype output was read.

### Copy number: one deviation in application, not in the rule

The rule's formula assumed every target is copy-neutral. That held on the new-diagnosis samples but not here:

| locus class | targets | treatment |
|---|---|---|
| chrX: AML5 is XY (DRAGEN ploidy); AML4's ploidy metrics are empty, but its chrX VAFs (0.68, 0.79) exceed p/2, so the locus is hemizygous (male) or LOH (female) | AML4 T1, T6; AML5 T16 | CCF = VAF / p. **Same verdict** as 2·VAF/p. |
| inside a clonal single-copy loss (DRAGEN CNV PASS, CN 1, segment mean 0.56, which matches CN 1 at p = 0.88) | AML5 T4, T6 | CCF = VAF·(p + 2(1−p)) / p, the standard CCF with total CN 1 and multiplicity 1. T6 has the **same verdict** either way. **T4 differs**: 2·VAF/p = 1.07, copy-number-aware = 0.59 → **subclonal**. |
| VAF > p/2 with no CNV call (probable copy-neutral LOH) | AML5 T1, sAML1 T11 | Clonal under 2·VAF/p and under any plausible copy state: VAF/p ≥ 0.8 already implies CCF ≥ 0.8 for hemizygous, LOH or single-copy-gain states. |

2·VAF/p is the copy-neutral special case of the copy-number-aware formula. Using the copy-aware value where the copy state is known is therefore the rule applied faithfully, not a new rule. **AML5 T4 is the only target where the choice changes the outcome**, and it is recorded as subclonal.

## Result

| Sample | p | targets | Clonal (CCF) | Subclonal (CCF) | Not PASS / absent |
|---|---|---|---|---|---|
| AML4 | 0.77 | 11 | **T1 (0.88), T2 (1.04), T3 (0.88), T4 (0.85), T6 (1.03), T7 (1.00), T8 (1.14), T9 (0.96), T10 (1.10)** | T5 0.57, T11 0.49 | – |
| AML5 | 0.88 | 16 | **T1 (LOH, ≥ 0.9), T2 (1.00), T3 (1.01), T5 (0.94), T6 (0.97), T7 (1.11), T9 (0.88), T10 (0.99), T12 (1.08)** | T4 0.59 (CN 1), T8 0.74, T11 0.75, T16 0.18 | T13, T15 (no call), T14 (filtered: PoN, MQ) |
| sAML1 | 0.67 | 16 | **T2 (1.11), T3 (1.06), T4 (1.12), T5 (1.00), T7 (0.99), T9 (1.10), T10 (0.97), T11 (LOH, ≥ 0.9), T12 (0.84), T13 (1.00), T14 (0.96)** | T1 0.51, T6 0.56, T8 0.53, T15 0.47, T16 0.34 | – |

Clonal sets: **AML4: T1–T4, T6–T10 (9). AML5: T1–T3, T5–T7, T9, T10, T12 (9). sAML1: T2–T5, T7, T9–T14 (11).**
No two clonal targets in a sample lie within the same molecule (all are in different genes), so each counts as its own locus.

Purity sanity check: the genome-wide PASS somatic VAF histogram (tumour DP ≥ 50) has its upper mode at 0.3–0.5 (AML4, AML5) and 0.3–0.4 (sAML1). That agrees with p = 0.77 / 0.88 / 0.67.

## Margins: narrower than on the new-diagnosis samples

WES depth at targets is 23–584× (median ≈ 50× in AML4/AML5, ≈ 250× in sAML1), so the 95 % Wilson bounds on CCF are wide:
- closest **subclonal** calls: AML5 T8 0.74 [0.49, 1.05] and T11 0.75 [0.50, 1.06];
- closest **clonal** calls: sAML1 T12 0.84 [0.65, 1.05], AML4 T4 0.85 [0.63, 1.10], AML4 T1 / T3 0.88, AML5 T9 0.88 [0.56, 1.25].

The rule is a point estimate and is applied as fixed. The pre-registration should still know that roughly six targets could flip under resampling. One option it may take: a sensitivity analysis that drops targets whose CCF bound straddles 0.8. That is the pass bar's call, not this ticket's.

## Consequence for the pool

AML4, AML5 and sAML1 move from calibration-only to **tumour-side evaluable**: six patients now carry a clonal set (PALM08E, 11E, 14E, AML4, AML5, sAML1). PALM07E and PALM21E stay calibration-only.

## Reproduction

HPC working folder `q87_clonal/`: `wes_old.py` (sbatch job 9066038, `short`, 2 GB, under 1 min) dumps the per-target tumour/normal AD/AF, the CNV overlaps and the genome-wide PASS VAF histogram. `ccf.py` (log `ccf.log` there) computes p, CCF and bounds. The README section "Clonal targets, old-panel restart samples (ticket 87)" holds the WES paths.
