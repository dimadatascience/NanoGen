# PALM09E indel MUT loss in v3 consensus counts

> Wayfinder ticket *PALM09E indel MUT loss in v3 consensus counts* (#70), map *T/B negative controls vs blast calls* (#60). It feeds the verdict ticket (#69).

**Patient data:** aggregate numbers only, targets as T1…Tn. The scripts and the full tables are on the HPC, in the working folder under `q70_indel/` (`indel_reads.py`, `indel_reads2.py`, `out/`). Slurm jobs 9002992 and 9002999 ran on `short`, with pysam 0.22.1 in the `acox1/scmseq:1.0` image.

## Question

At PALM09E T1, a 4 bp insertion, v3 calls MUT in 34 of 1,428 covered Leukemic cells; `main` calls 1,272. At T2, a 1 bp insertion, the counts are 10 vs 72. Is this:

- (a) indel representation in PALM09E's BAM;
- (b) labelling in v3's CONSENSUS, where `main` and v3 label the same molecules differently; or
- (c) a real absence?

The comparator is PALM22E T4, the same insertion as PALM09E T1, where v3 separates the tumour cleanly.

## Answer in brief

**Mostly (b) with a weak, genuine signal underneath. Not (a).** v3 is not mislabelling the molecules: PALM09E's reads do not carry a tumour-specific insertion at a level that can be called per UMI.

1. **Representation is not the cause.** Both BAMs come from minimap2 2.28. At T1, 88 % of length-compatible insertions sit exactly on the BED anchor (PALM22E T4: 94.5 %). The shifted remainder is identical in Leukemic and negative cells. Matching tolerant of position (any offset within ±12 bp, length ±1) moves the Leukemic cells with a MUT-like UMI from 52 to 92, and negatives from 6 to 9. That falls far short of `main`'s 1,272.
2. **The `main`-vs-v3 gap is the call rule, not the molecules.** At an indel of > 2 bp, `main` labels a UMI MUT when **more than 3 reads** carry an alt-like entry, whatever the fraction (substring match). v3 needs ≥ 3 ALT reads **and** ≥ 60 % of the usable reads. At PALM09E T1, 16 % of reads in *negative* cells carry the exact insertion. Under `main`'s rule, 60 % of Leukemic UMIs and 53 % of negative UMIs are MUT, and `main` calls 637 of 956 covered negative cells. `main`'s 1,272 Leukemic calls are this background, not the clone.
3. **The tumour signal is real but weak, and it sits in mixed UMIs.** The excess of exact-insertion reads in Leukemic over negative cells is **+4.1 points** at PALM09E T1 (+5.8 with MAPQ < 20 reads included), against **+14.4 points** at PALM22E T4. Only 1.4 % of PALM09E's Leukemic UMIs reach ≥ 60 % insertion reads (negatives 0.6 %), against 13 % at PALM22E T4 (negatives 2 %). No per-UMI majority rule can turn this into per-cell calls.
4. **T2 is (c) in these data:** there is no Leukemic-vs-negative difference in any representation (exact-insertion reads 9.3 % vs 9.1 %; identical UMI and offset distributions).

**Consequence for the verdict:** PALM09E has no target with a usable tumour MUT signal, on either count source. PALM09E T1 does not fail because of a v3 defect. It cannot be read as a clonal target at all, and its criteria 2 and 3a outcomes are uninformative about v3. The verdict ticket decides how to record that (see *Implications*).

## Method

Cells: QC cells grouped as in the full-run coverage (#66): Leukemic, Ambiguous, negatives = T/B/NK. Two analyses:

- **A. Count cross-tab.** Per-cell WT/MUT UMIs and calls at T1/T2 from the `main` and v3 full-run outputs (#66).
- **B. Read level.** Primary reads in the canonical tagged BAM (#62) that cover the BED anchor. Each read with MAPQ ≥ 20 (both pipelines' filter) is classed by the insertion nearest the anchor:
  - **exact (v3's rule):** at the anchor, length within ±1 bp, prefix agreeing;
  - **anchor, sequence mismatch;**
  - **shifted:** length-compatible insertion at another offset within ±12 bp;
  - another insertion, a deletion only, a soft clip near the anchor, or no indel.

  UMIs are raw `UB` values, without v3's Levenshtein merge, so these UMI counts are close to v3's but not identical (52 vs 34 Leukemic cells with a MUT UMI at T1 under v3's rule). A follow-up job tabulated the exact-insertion fraction by MAPQ band and the per-UMI distribution of the exact-insertion read fraction.

## Results

### A. `main` vs v3 counts and calls (QC cells)

| target | group | covered cells (v3 / main) | v3 WT / MUT UMIs | main WT / MUT UMIs | MUT calls v3 / main | main-only MUT |
|---|---|---|---|---|---|---|
| PALM09E T1 | Leukemic | 1,428 / 1,447 | 8,017 / 34 | 11,884 / 4,195 | 34 / 1,272 | 1,240 |
| PALM09E T1 | negatives | 956 / 970 | 3,301 / 3 | 5,051 / 1,224 | 3 / 637 | 634 |
| PALM09E T2 | Leukemic | 785 / 988 | 2,906 / 10 | 19,228 / 87 | 10 / 72 | 71 |
| PALM09E T2 | negatives | 678 / 885 | 2,355 / 8 | 14,280 / 60 | 8 / 54 | 53 |
| PALM22E T4 | Leukemic | 4,237 / 4,460 | 22,887 / 4,738 | 33,765 / 40,491 | 2,514 / 4,386 | 1,877 |
| PALM22E T4 | negatives | 2,029 / 2,055 | 9,491 / 57 | 8,311 / 7,965 | 55 / 1,957 | 1,902 |

`main`'s MUT calls are near-universal in negatives at both 4 bp-insertion targets, so its rule does not discriminate there. v3's calls are almost all a subset of `main`'s (2 / 5 v3-only cells).

### B1. Where the insertion sits (reads covering the anchor, MAPQ ≥ 20)

| target | group | reads | MAPQ < 20 (dropped) | exact | anchor, seq mismatch | shifted | other ins | del only | no indel |
|---|---|---|---|---|---|---|---|---|---|
| PALM09E T1 | Leukemic | 321,860 | 34.6 % | **19.9 %** | 1.1 % | 2.7 % | 9.9 % | 10.4 % | 55.4 % |
| PALM09E T1 | negatives | 126,278 | 32.4 % | **15.9 %** | 0.9 % | 2.2 % | 9.6 % | 11.1 % | 59.7 % |
| PALM09E T2 | Leukemic | 2,208,527 | 0.0 % | 9.3 % | 0.6 % | 4.2 % | 1.0 % | 38.8 % | 46.0 % |
| PALM09E T2 | negatives | 1,525,182 | 0.1 % | 9.1 % | 0.6 % | 4.2 % | 1.0 % | 38.7 % | 46.4 % |
| PALM22E T4 | Leukemic | 2,449,657 | 7.7 % | **38.0 %** | 0.7 % | 2.1 % | 7.1 % | 6.4 % | 45.7 % |
| PALM22E T4 | negatives | 508,896 | 6.2 % | **23.6 %** | 0.5 % | 1.6 % | 6.9 % | 8.2 % | 59.3 % |

Length-compatible insertions at the anchor offset (0) vs elsewhere:

- PALM09E T1: 0.879 Leukemic, 0.878 negatives; the largest other offset is −2, at 0.055.
- PALM22E T4: 0.945 Leukemic, 0.935 negatives.
- PALM09E T2 (homopolymer): 0.655 Leukemic, 0.651 negatives; the remainder spreads over ±12 bp in both groups alike.

### B2. MAPQ: a BAM-level difference that does not create the gap

| target | group | MAPQ 1–19: reads, exact-insertion frac | MAPQ ≥ 20: reads, exact-insertion frac |
|---|---|---|---|
| PALM09E T1 | Leukemic | 164,464, 0.614 | 321,860, 0.200 |
| PALM09E T1 | negatives | 58,306, 0.549 | 126,278, 0.159 |
| PALM22E T4 | Leukemic | 202,696, 0.492 | 2,449,657, 0.380 |
| PALM22E T4 | negatives | 33,122, 0.231 | 508,896, 0.236 |

Insertion-bearing reads are pushed to low MAPQ in PALM09E far more than in PALM22E: 34 % of reads at the anchor are MAPQ < 20, against 8 %. This happens in negatives as much as in Leukemic cells, and both pipelines drop these reads. The Leukemic excess inside the low-MAPQ band is also small (+6.5 points). Over all MAPQ, the excess is +5.8 points against PALM22E T4's +15.

### B3. Per raw UMI (≥ 3 reads, MAPQ ≥ 20): exact-insertion read fraction

| target | group | UMIs | 0 | (0, 0.2) | [0.2, 0.4) | [0.4, 0.6) | [0.6, 0.8) | [0.8, 1] |
|---|---|---|---|---|---|---|---|---|
| PALM09E T1 | Leukemic | 9,891 | 0.092 | 0.456 | 0.398 | 0.039 | 0.012 | 0.002 |
| PALM09E T1 | negatives | 3,716 | 0.115 | 0.612 | 0.252 | 0.015 | 0.006 | 0.000 |
| PALM09E T2 | Leukemic | 19,854 | 0.526 | 0.256 | 0.191 | 0.019 | 0.008 | 0.000 |
| PALM09E T2 | negatives | 14,092 | 0.534 | 0.247 | 0.189 | 0.020 | 0.009 | 0.000 |
| PALM22E T4 | Leukemic | 55,926 | 0.052 | 0.111 | 0.363 | 0.346 | 0.111 | 0.018 |
| PALM22E T4 | negatives | 11,639 | 0.086 | 0.227 | 0.633 | 0.035 | 0.017 | 0.003 |

In negative cells, a typical UMI at a 4 bp-insertion target carries 20–40 % insertion reads, on both samples. A single molecule cannot do that. This is consistent with intermolecular mixing after UMI tagging (PCR chimeras / template switching), which spreads mutant sequence across families. The background explains `main`'s calls in negatives. At PALM22E T4, the Leukemic distribution moves to [0.4, 0.8), and v3's ≥ 60 % rule catches the upper part. At PALM09E T1, it moves only from (0, 0.2) to [0.2, 0.4), below any majority rule.

### Per cell (cells with ≥ 1 UMI of ≥ 3 reads): cells with a MUT-like UMI

| target | group | cells | exact rule (v3) | position-tolerant |
|---|---|---|---|---|
| PALM09E T1 | Leukemic | 1,448 | 52 | 92 |
| PALM09E T1 | negatives | 968 | 6 | 9 |
| PALM09E T2 | Leukemic | 855 | 32 | 138 |
| PALM09E T2 | negatives | 800 | 27 | 93 |
| PALM22E T4 | Leukemic | 4,455 | 2,739 | 3,273 |
| PALM22E T4 | negatives | 2,049 | 89 | 104 |

## Interpretation

- **(a) Representation: no.** Insertion placement matches PALM22E and does not differ by group. The one real BAM difference, many more low-MAPQ reads at T1, affects both groups alike and carries little tumour excess.
- **(b) Labelling: the rules differ, but v3 is the one that reads the data right.** `main`'s "> 3 alt reads" rule labels the mixing background as MUT, which explains its 1,272 vs 34. v3's majority rule is doing its job. It would fail only if the clone's mutant molecules had a within-UMI insertion fraction below 0.6, and PALM22E T4 shows they do not: the same allele there gives 13 % of tumour UMIs at ≥ 0.6.
- **(c) Real absence: partly.** At T1, a weak tumour excess exists: a quarter to two-fifths of PALM22E's in reads, and under a tenth in excess clean MUT UMIs. It fits a minority of Leukemic cells carrying T1, or low mutant expression in this sample. The data cannot tell these apart, and the earlier π̂ attempts (#63) failed for the same reason. At T2 there is no excess.

A side observation, outside this ticket's question: at PALM22E T4, position-tolerant matching would add ~19 % more Leukemic cells with a MUT UMI (2,739 → 3,273) at little cost in negatives (89 → 104). That bears on the fog item *Alternative filter if v3 fails*. It is not a change made here (code changes are out of scope for the map).

## Implications for the verdict (#69)

- PALM09E T1 should not be read as a clonal target for criteria 2 and 3a. Neither count source yields a tumour signal there, and the failure is not attributable to v3. Under the pre-registration's own logic, this makes PALM09E's clonal criteria **not evaluable** rather than *failed*. Whether that re-labelling is allowed after pre-registration is the verdict ticket's call.
- PALM09E's criterion 1 (calibration) is unaffected and stands as scored.
- `main`'s indel MUT calls at 4 bp-insertion targets are not a valid sensitivity benchmark: they are near-universal in negatives.
