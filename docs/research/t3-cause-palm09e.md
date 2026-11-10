# PALM09E T3 cause: germline, artefact or B-cell over-dispersion

Wayfinder ticket: *PALM09E T3 cause: germline, artefact or B-cell over-dispersion* (#73), map *Confirm genotype-v3 vs v3-ub on new diagnosis samples* (#72).

- **Pre-registration:** [`t3-cause-preregistration.md`](t3-cause-preregistration.md), committed at `8b8e330`, with its A4 threshold amended at `ba81281` before any data look.
- **Code:** [`t3-cause/t3_cause.py`](t3-cause/t3_cause.py), run on `hpc_new` (job 9019141, `short`), in the pinned NanoGen image (numpy 2.0.1, scipy 1.14.0, pysam 0.22.1, pandas 2.2.2).
- **Outputs:** aggregate only, in [`t3-cause/out/`](t3-cause/out/).
- **Inputs:** the genotype-v3 full run of PALM09E (#66) and the canonical tagged BAM (#62).

## Verdict: inconclusive

The pre-registered rule returns **inconclusive** (§4, row 5): both an artefact signal (A3) and an over-dispersion signal (D1 and D2) hold. Under §5 this means:

- **v3 and v3-ub both go to the new samples unchanged.** T3 neither supports nor rules out v3-ub.
- **The PALM09R pair check graduates from the fog.**

## Checks

- The self-test passed. On a known over-dispersed group, D1 gives p = 4 × 10⁻⁹; on a binomial group it gives p = 0.49.
- The Neg pooled rate recomputed from the counts equals the run's c_s at T3 (0.0396).
- Neg is exactly T ∪ NK ∪ B, and the three sets don't overlap.

## Signals

| signal | measurement | threshold | holds |
|---|---|---|---|
| **G** germline | UMI MUT rate (Wilson 95 % CI): T 3.4 % (2.3–5.1), B 13.5 % (5.9–28.0), Leuk 13.0 % (10.3–16.4) | every CI overlaps [0.3, 0.7], or every lower bound ≥ 0.9 | no |
| **A1** strand | Neg reads: f₊(ALT) 0.54 vs f₊(REF) 0.36, Δ = 0.19, p < 10⁻³⁰⁰ | p < 0.01 and \|Δ\| ≥ 0.25 | no (Δ below the effect bar) |
| **A2** read position | Neg reads: median distance to the nearer read end is 188 (ALT) vs 94 (REF), ratio 2.0 | p < 0.01 and ratio ≤ 0.5 | no (ALT sits *farther* from the ends) |
| **A3** molecule purity | median ALT fraction of ALT-majority raw molecules: Neg 0.67 (n = 290) vs Leuk 1.00 (n = 2,153); difference 0.33, p = 2 × 10⁻¹¹³ | p < 0.01 and difference ≥ 0.15 | **yes** |
| **A4** context | Neg raw ALT rate at T3 is 6.0 %; same-substitution background over 33 positions has median 0.26 % and 95th percentile 0.76 % (ratio 8.0) | ratio ≤ 1 | no (T3 is an outlier) |
| **D1** over-dispersion | Neg cells: beta-binomial ρ̂ = 0.22, LR = 12.9, p = 1.7 × 10⁻⁴ | p < 0.05 | **yes** |
| **D2** B-specificity | covered cells with ≥ 1 MUT UMI: B 5/20 vs T/NK 26/292, one-sided Fisher p = 0.037 | p < 0.05 | **yes** |
| **U** homogeneity | T/NK, B, Leuk: χ² = 41.2, Monte Carlo p = 2 × 10⁻⁵ | p ≥ 0.05 | no |

**Germline is ruled out** under any reading: every negative lineage reads MUT at a few per cent, not around 50 %, and the rate is far from uniform across lineages.

## What the data show (reported, not deciding)

These observations don't change the verdict. They set out what the PALM09R check has to separate.

**1. Two distinct populations of ALT signal in negatives.**

- **In raw reads, T3 is a strong site-specific outlier.** Negatives read ALT in 6.0 % of 2.8 M reads, against a same-substitution background whose maximum is 1.3 %.
- **Almost none of that signal forms whole mutant molecules.** Of 11,184 raw negative molecules (≥ 3 reads):
  - 2,573 carry ALT in under 20 % of their reads;
  - 1,716 carry it in 20–50 %;
  - 262 in 50–80 %;
  - only 28 (0.25 %) in ≥ 80 %.
- **Tumour cells look different.** 1,798 of 6,845 tumour molecules (26 %) are ≥ 80 % ALT, and the median ALT-majority molecule is pure.
- **A vagrant mutant transcript should arrive as a whole mutant molecule**, as it does in tumour cells. Partial-ALT molecules are what a site-specific read-level error produces, or UMI collision with tumour molecules.
- **Against a plain basecalling error:**
  - ALT bases have *higher* base quality than REF bases (median 31 vs 27);
  - ALT reads sit farther from read ends;
  - ALT reads are moderately enriched on the + strand.

  So the negative-cell ALT reads come from a read population that differs from that of the REF reads. In tumour cells, ALT and REF reads are alike (distance 52 vs 53, Δ strand 0.03).

**2. The over-dispersion is not "a few cells".**

- 31 negatives carry MUT UMIs: 30 have exactly one, and one has two.
- D1's excess therefore isn't concentration of MUT in a few cells. It is consistent with the one-UMI pattern already reported in #67 (MUT in 8.0 % of one-UMI negatives vs a pooled 4.0 %): MUT UMIs turn up more often in low-depth cells than a constant per-UMI rate predicts. This was not tested here.
- MUT carriers are not doublets (0 % flagged in both carriers and others), and they are no deeper than other negatives (median total UMIs 6 vs 7, p = 0.79).

**3. B-specificity is weak and target-local.**

- At T3, B cells read MUT 3.9× as often as T/NK, at p = 0.037 on 20 cells.
- At the two other targets with ≥ 20 covered B cells (T1, T2), B cells read **no** MUT UMIs (0/83 and 0/140), while T/NK read 0.09 % and 0.36 %.
- Vagrant contamination that favours B cells would show at those targets too. T4 has too few B cells to check.

**4. Editing caveat.** The target gene is on the minus strand (public annotation; no GTF was available on the HPC), so the genomic G>A reads as C>U on the transcript, an APOBEC-type change. RNA data can't separate C>U editing from a DNA change. Editing would produce whole edited molecules, though, which the negatives mostly lack.

**Overall reading (not a verdict).** The evidence leans towards a **site-specific read-level artefact** at T3 in negative cells, plus a small residual of whole mutant molecules. That is closer to "germline or artefact" (outlier flag under-sensitive) than to B-cell vagrant over-dispersion. The leaning matters for the PALM09R check: an artefact of this kind should recur at relapse with the same partial-molecule profile.

## Deviations

1. **A4 threshold amended before any data look** (`ba81281`), as disclosed in the pre-registration.
2. **The trinucleotide-matched A4 background is empty.** No position within ±100 bp shares T3's trinucleotide context and has ≥ 200 negative reads, so only the all-G background is reported. A4 doesn't hold on it either way.
3. **The transcript strand comes from public annotation**, not from a GTF on the HPC (none found under the resources folder).
