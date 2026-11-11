# Regenerated census tagged BAMs: short-read known barcode list

Wayfinder ticket #85 (map #72), after the adopt decision in #82
([`known-bc-list-palm14e.md`](known-bc-list-palm14e.md)).

**Outcome:**
- The regenerated BAM becomes canonical for PALM08E, PALM11E, PALM14E and PALM21E.
- PALM07E keeps its canonical library-1 BAM, because its listed FASTQ is a library-2 enrichment.

Aggregate counts only. Targets are T1…Tn in BED order. Paths, barcodes and the target key stay in the HPC
README (working folder `q85_regen/`).

## Setup

- **Code:** prototype `aa30310` (`prototype/known-bc-list`), pinned image `sha256:764368d6…` (`blaze2` 2.4.0).
  - A driver (`regen.nf`, HPC only) sits in a copy of the prototype tree and imports its own modules:
    BLAZE `--no-whitelisting --known-bc-list` → MERGE_READS → MINIMAP2 → FIXTAGS, then `samtools sort` and
    `index`.
  - It stops before SPLIT_BAM: the per-cell tail runs later from the tagged BAM.
- **Known list:** the sample's short-read QC list; for PALM07E, the derived library-1 list.
- **Lanes:** every FASTQ in the sample-paths sheet, one lane per file.

  | Sample | Lanes | Reads |
  |---|---|---|
  | PALM07E | 1 | 38.2 M |
  | PALM08E | 1 | 90.7 M |
  | PALM11E | 1 | 30.8 M |
  | PALM21E | 2 | 35.5 M + 26.6 M |

- **PALM14E:** the #82 pilot's known-list BAM (both lanes), sorted and indexed the same way.
- **Coverage, outcome-blind:** the power check's BAM-level proxy (`power_check.py coverage` + `summarise`),
  unchanged and with the same parameters. It reads no allele.
  - It ran on each regenerated BAM. The canonical BAM's numbers are the power check's own tables, from the
    same code and image.
  - `neg` = covered T/B/NK negatives in QC; `tum` = covered Leukemic cells.
- **Gate (posted on the ticket before any regenerated BAM existed):** pooled over the intake's well-covered
  targets, regenerated `neg` ≥ canonical and regenerated `tum` ≥ canonical.
  - Well-covered targets: all targets, except PALM08E T8–T11, PALM11E T2 and PALM21E T6.

## Result

| Sample | neg canonical | neg regenerated | tum canonical | tum regenerated | Canonical BAM |
|---|---|---|---|---|---|
| PALM07E | 633 | 84 | 1,233 | 63 | **unchanged** (library 1) |
| PALM08E | 361 | 431 | 7,855 | 7,932 | regenerated |
| PALM11E | 282 | 366 | 5,766 | 6,841 | regenerated |
| PALM14E | 194 | 272 | 1,890 | 2,256 | regenerated |
| PALM21E | 1,341 | 1,857 | 3,426 | 3,536 | regenerated |

The excluded targets have zero coverage in both BAMs, so pooling over all targets gives the same sums. In the
four regenerated samples, no target falls below the canonical BAM's `neg` or `tum`.

### Per target

| Sample | Target | neg canon | neg regen | tum canon | tum regen |
|---|---|---|---|---|---|
| PALM08E | T1 | 2 | 6 | 16 | 16 |
| PALM08E | T2 | 0 | 0 | 30 | 30 |
| PALM08E | T3 | 86 | 94 | 2,447 | 2,463 |
| PALM08E | T4 | 96 | 121 | 1,934 | 1,958 |
| PALM08E | T5 | 36 | 43 | 855 | 860 |
| PALM08E | T6 | 1 | 1 | 88 | 91 |
| PALM08E | T7 | 140 | 166 | 2,485 | 2,514 |
| PALM11E | T1 | 156 | 231 | 3,310 | 4,179 |
| PALM11E | T3 | 13 | 14 | 1,065 | 1,084 |
| PALM11E | T4 | 20 | 21 | 214 | 253 |
| PALM11E | T5 | 93 | 100 | 1,177 | 1,325 |
| PALM14E | T1 | 10 | 16 | 92 | 135 |
| PALM14E | T2 | 5 | 16 | 76 | 117 |
| PALM14E | T3 | 3 | 7 | 108 | 131 |
| PALM14E | T4 | 2 | 4 | 85 | 108 |
| PALM14E | T5 | 2 | 5 | 87 | 109 |
| PALM14E | T6 | 87 | 120 | 762 | 875 |
| PALM14E | T7 | 85 | 104 | 680 | 781 |
| PALM21E | T1 | 89 | 132 | 255 | 264 |
| PALM21E | T2 | 78 | 98 | 202 | 218 |
| PALM21E | T3 | 3 | 3 | 7 | 7 |
| PALM21E | T4 | 19 | 41 | 39 | 42 |
| PALM21E | T5 | 19 | 24 | 95 | 100 |
| PALM21E | T7 | 182 | 199 | 267 | 272 |
| PALM21E | T8 | 267 | 317 | 759 | 775 |
| PALM21E | T9 | 106 | 158 | 312 | 320 |
| PALM21E | T10 | 152 | 241 | 526 | 538 |
| PALM21E | T11 | 326 | 507 | 771 | 797 |
| PALM21E | T12 | 100 | 137 | 193 | 203 |

## Findings

- **PALM07E's listed FASTQ is library 2.** It comes from a nanopore run sequenced after the canonical library-1
  BAM was made. Evidence:
  - **Mock** (first 250k reads): 1,135 tagged mapped reads with the library-1 list, against 109,287 with the
    original QC list. The original QC list is library 2 (#74), and 109,287 is in line with the other samples'
    mocks (112–152 k).
  - **Full run** with the library-1 list: 84 covered negatives and 63 tumour cells, against 633 and 1,233 for
    the canonical BAM.
  - So PALM07E stays on library 1. Whether a library-2 PALM07E BAM would serve better is a separate question;
    it needs library-2 labels and negatives.
- **The redemux / reseq concern didn't materialise.** PALM11E (canonical from a separate redemux) gains 30 %
  negatives and 19 % tumour cells from the listed FASTQ. PALM21E (two FASTQs) gains 38 % negatives and 3 %
  tumour cells.
- **The gains are mostly in negatives.** That's the scarce input, consistent with the #82 pilot. Tumour gains
  are small except in PALM11E and PALM14E.

## Resources and deviations

- **Mock first:** 250k reads per lane through the whole driver and the proxy. No failed task.
- **Full runs:** every head exited 0.
- **Deviation:** the sourced per-sample resource env didn't reach the Nextflow heads (`SBATCH_EXPORT=NONE` in
  the submitting shell). Tasks therefore requested `regen.config`'s defaults:

  | Process | Requested (default) | Pre-computed request |
  |---|---|---|
  | BLAZE | 8 GB / 60 min | 7 GB |
  | MINIMAP2 | 34 GB / 60 min | 34 GB |
  | SORT | 14 GB / 60 min | 12 GB |
  | MERGE | 2 GB | — |

- **PALM08E retries:** BLAZE hit 60 min and completed on retry (16 GB / 2 h requested; 1 h 26 m, 6.4 GB peak).
  FIXTAGS hit 20 min and completed on retry in 21 m 41 s. The retry's 16 GB exceeds measured + 50 %.
  - For a future run: set `--export=ALL`. Also, BLAZE scales with reads, not file size: PALM08E's peak is
    6.4 GB, against 3.1–3.9 GB elsewhere.
- **Measured peaks (requests for the full runs):**

  | Process | Peak | Wall |
  |---|---|---|
  | BLAZE | 3.1–6.4 GB | 17–86 min |
  | MINIMAP2 | 15–22.6 GB | 17–54 min |
  | SORT | 8.3–8.4 GB | ≤ 6 min |
  | Proxy | ≤ 1.2 GB | ≤ 12 min |
