# Short-read known barcode list vs BLAZE whitelisting: PALM14E pilot

Wayfinder ticket #82 (map #72). Decision: **adopt**. Regenerate the census tagged BAMs from FASTQ with the
short-read QC list as BLAZE's known barcode list.

Aggregate counts only. Targets are T1…T7 in BED order; the key, paths and barcodes stay on the HPC
(working folder `q82_knownbc/`).

## Question

BLAZE keeps a barcode only if its high-quality read count is at least `quantile(top expect_cells, 0.95) / 20`,
set per lane. On a targeted panel this may drop cells with few reads at the targets. Does assigning reads
to the short-read cell barcodes (`--no-whitelisting --known-bc-list`) give enough covered cells per target to
replace the current barcode step?

## Setup

- **Code:** `prototype/known-bc-list` @ `aa30310` (genotype-v3 `61d8809` + an optional `known_bc` samplesheet
  column). Pinned image `sha256:764368d6…`, `blaze2` 2.4.0. One SHA and one image for all arms.
- **Sample:** PALM14E, both nanopore FASTQs (barcode28, barcode29) as two lanes.
- **Arms:**
  - **default:** BLAZE `--expect-cells 5000 --max-edit-distance 1`;
  - **known:** `known_bc` = the PALM14E short-read QC list (5,945 cells);
  - **canonical:** the current canonical tagged BAM, through the tagged-BAM entry point.
- **Known list = QC list, not Cell Ranger filtered barcodes:** no filtered lists sat next to the census inputs,
  and every scored cell is in the QC list.
- **Metrics (outcome-blind):** QC cells only; covered = WT + MUT ≥ 1 UMI. Only the WT + MUT sum was read: no MUT
  count, MUT fraction or call.
- **Adoption rule:** pre-registered on the ticket before any run.
- **Runs:** heads completed with exit 0. Every task ran on its first attempt, with no failure or retry.
  Wall time: default 3 h 35 m, known 2 h 04 m, canonical 1 h 37 m.

## Coverage per target

Covered QC cells by group (neg = T/B/NK negatives, tum = Leukemic, amb = Ambiguous) and WT + MUT UMIs.

| Target | neg default | neg known | neg canonical | tum default | tum known | tum canonical | amb default | amb known | amb canonical | UMIs default | UMIs known | UMIs canonical |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T1 | 10 | 14 | 8 | 92 | 99 | 73 | 33 | 36 | 31 | 304 | 319 | 235 |
| T2 | 11 | 15 | 5 | 116 | 117 | 76 | 19 | 22 | 16 | 191 | 205 | 128 |
| T3 | 2 | 4 | 2 | 94 | 97 | 86 | 15 | 17 | 15 | 357 | 364 | 283 |
| T4 | 2 | 4 | 2 | 93 | 95 | 81 | 13 | 15 | 11 | 263 | 269 | 204 |
| T5 | 2 | 4 | 2 | 92 | 94 | 73 | 13 | 15 | 13 | 251 | 258 | 197 |
| T6 | 89 | 116 | 84 | 774 | 825 | 718 | 119 | 133 | 107 | 2,911 | 3,014 | 2,191 |
| T7 | 73 | 82 | 63 | 602 | 636 | 553 | 96 | 105 | 87 | 1,695 | 1,746 | 1,402 |
| **pooled** | **189** | **239** | **166** | **1,863** | **1,963** | **1,660** | 308 | 343 | 280 | 5,972 | 6,175 | 4,640 |

Known / default ratios:

| Target | neg | tum | UMIs |
|---|---|---|---|
| T1 | 1.400 | 1.076 | 1.049 |
| T2 | 1.364 | 1.009 | 1.073 |
| T3 | 2.000 | 1.032 | 1.020 |
| T4 | 2.000 | 1.022 | 1.023 |
| T5 | 2.000 | 1.022 | 1.028 |
| T6 | 1.303 | 1.066 | 1.035 |
| T7 | 1.123 | 1.056 | 1.030 |
| pooled | **1.265** | **1.054** | **1.034** |

Cells in the per-cell output: 11,305 in default (4,100 of them QC); 5,945 in known (all QC); 4,997 in
canonical (4,095 QC).

## Adoption rule (pre-registered) applied

| Criterion | Threshold | Observed | Pass |
|---|---|---|---|
| (a) pooled `neg` known/default | ≥ 1.10 | 1.265 | yes |
| (b) pooled `tum` known/default | ≥ 0.98 | 1.054 | yes |
| (b) pooled `umi` known/default | ≥ 0.98 | 1.034 | yes |
| (c) every target, `neg` and `tum` | ≥ 0.95 | min 1.123 (neg T7), 1.009 (tum T2) | yes |
| (d) known ≥ canonical, pooled `neg` and `tum` | ≥ | 239 ≥ 166; 1,963 ≥ 1,660 | yes |

**Adopt.**

## Sanity

| Arm | On-target reads, QC barcodes | On-target reads, non-QC barcodes | QC CBs with on-target reads | Non-QC CBs |
|---|---|---|---|---|
| default | 6,444,566 | 3,204,757 | 4,071 | 2,830 |
| known | 6,478,198 | 0 | 5,768 | 0 |
| canonical | 3,811,311 | 1,905,913 | 4,065 | 843 |

- **Mis-assignment is bounded.** Reads from non-QC barcodes could be pulled into a QC cell at edit distance 1. The
  known arm has only 33.6 k (0.5 %) more on-target reads in QC barcodes than default, while covered negatives rise
  26 %. The gain is cells with few reads that the cutoff dropped, not reads moved from neighbours.
- **The per-lane cutoff is unstable.** BLAZE's default whitelist keeps 4,997 barcodes in lane 28 (4,095 in the QC
  list) and 9,750 in lane 29 (only 2,754 in the QC list).
- **The canonical BAM is lane 28 alone.** It has the same 4,997 barcodes as lane 28's default whitelist. Part of
  canonical's deficit against default (pooled UMIs 0.78×) is the missing second lane. That part isn't the barcode
  step.

## Caveats

- **One pilot sample.** The negative gain on T1–T5 rests on 2–15 cells per target. The pooled gain comes mostly
  from T6 and T7.
- **Low absolute negative coverage.** Even with the known list, PALM14E has no target with ≥ 200 covered negatives
  (the evaluability bar carried from the T/B negative-controls map). This is for the power check to judge.
- **Regeneration starts from the FASTQs in the sample-paths sheet.** Where a canonical BAM came from a different
  demultiplexing of the nanopore run (PALM11E `redemux`, PALM21E `reseq`), the FASTQ may not hold the same reads.
  Each regenerated BAM is checked against its canonical BAM before it replaces it.
