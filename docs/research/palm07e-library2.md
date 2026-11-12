# PALM07E: library 2 from the October FASTQ

Ticket: [PALM07E: switch to library 2 from the October FASTQ?](https://github.com/dimadatascience/NanoGen/issues/88), map [#72](https://github.com/dimadatascience/NanoGen/issues/72).
Aggregate numbers only; targets T1–T5. Inputs, scripts and paths live in the HPC working folder (`q88_lib2/`, README section "ticket 88").

## Question

PALM07E has two 10x libraries. The census scored library 1, because no library-2 BAM existed. Its listed FASTQ (a later nanopore run) turned out to be a library-2 enrichment ([#85](https://github.com/dimadatascience/NanoGen/issues/85)). Is a library-2 BAM regenerated from that FASTQ at least as well covered as library 1?

## Method (gate fixed in the ticket before the run)

- **Inputs.** From the full Seurat annotation, same rule as the library-1 intake ([#74](https://github.com/dimadatascience/NanoGen/issues/74)): QC universe = library-2 cells, negatives = T/NK ∪ B, labels = `ct`.
- **Regeneration.** The [#85](https://github.com/dimadatascience/NanoGen/issues/85) driver (prototype `aa30310`): BLAZE `--known-bc-list` → MERGE_READS → MINIMAP2 → FIXTAGS → sort + index, pinned image, one lane. Known list = the original QC list.
- **Coverage, outcome-blind.** The power check's BAM proxy (`power_check.py coverage` + `summarise`), unchanged. `neg` = covered negatives, `tum` = covered Leukemic cells.
- **Gate.** Pooled over T1–T5, library-2 `neg` ≥ 633 **and** `tum` ≥ 1,233 → switch to library 2.

## Input checks

| | Library 2 (derived) | Original file |
|---|---|---|
| QC cells | 6,637 | 6,637 (identical set) |
| Negatives (T/NK ∪ B) | 3,351 | 6,864 in the CSV; 3,351 of them in QC = the derived set exactly |
| Leukemic / Ambiguous | 2,352 / 934 | – |

The T/B/NK CSV is larger than the QC list because it also holds the library-1 negatives (3,513). Nothing else accounts for the difference.

## Result

| Target | neg lib1 | neg lib2 | tum lib1 | tum lib2 | E_B(c=1 %) lib1 | E_B(c=1 %) lib2 | evaluable lib1 → lib2 |
|---|---|---|---|---|---|---|---|
| T1 | 270 | 313 | 317 | 299 | 2.85 | 2.21 | yes → yes |
| T2 | 12 | 13 | 21 | 9 | 0.06 | 0.09 | no → no |
| T3 | 89 | 92 | 242 | 285 | 0.73 | 0.70 | no → no |
| T4 | 88 | 21 | 198 | 130 | 0.45 | 0.15 | no → no |
| T5 | 174 | 393 | 455 | 778 | 1.76 | 3.77 | no → **yes** |
| **Pooled** | **633** | **832** | **1,233** | **1,501** | **5.84** | **6.92** | |

**Gate passes** (832 ≥ 633, 1,501 ≥ 1,233). **PALM07E switches to library 2.**

- Pooled E_B at c = 1 % rises by 1.08 null calls (+18 %). T5 becomes evaluable.
- T4 loses most of its negatives (88 → 21) and T2 its tumour cells (21 → 9). Neither target was evaluable on library 1, and the gate is pooled. The final pool ([#86](https://github.com/dimadatascience/NanoGen/issues/86)) should use the library-2 per-target values.

## Runs and deviations

- **Resources** were sized from the library-1 run on the same FASTQ and from PALM11E scaled ×1.24, capped at ≤ measured + 50 %, with the env passed via `--export=ALL`.
- **MINIMAP2 was under-sized.** It hit its 32 min limit and completed on retry: 40 min, 22.6 GB peak.
- **Other steps:** BLAZE 29 min / 3.6 GB; FIXTAGS 6.6 min; SORT_INDEX 1.3 min / 8.3 GB; proxy 2 min.
- **BAM size:** 8.9 GB, vs 96 MB for the library-1 regeneration from the same FASTQ.
- **Full runs:** must pass the derived library-2 negatives CSV. v3 matches negatives by exact CB with no QC filter.
