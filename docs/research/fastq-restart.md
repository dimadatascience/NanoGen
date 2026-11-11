# FASTQ restart: power check on the old-panel FASTQ-only samples

Ticket: [FASTQ restart: power check on the old-panel FASTQ-only samples](https://github.com/dimadatascience/NanoGen/issues/84), map [confirm genotype-v3 vs v3-ub on new diagnosis samples](https://github.com/dimadatascience/NanoGen/issues/72).
Aggregate numbers only. Targets are T1…Tn in BED order, and the key is on the HPC (`q84_restart/`, README). No barcodes, paths, genes or coordinates appear here.

## Answer

AML4, AML5 and sAML1 **pass intake** and add **E_B = 9.3 expected null calls at c = 1 %** to the pool: AML4 2.5, AML5 2.3, sAML1 4.5. That is +41 % over the five-patient census (22.4).

**No old-panel target is evaluable on its own.** sAML1 T9 is borderline: 12/20 seeds raw, 0 at worst case. No target reaches E_B ≥ 5.7 at c = 1 %, the site screen's 80 % power point for a 3× miss.

A preview of the eight-patient pool (census on its canonical BAMs, plus the restart) gives:
- E_B 31.7 at c = 1 %, just past E\* = 30.9 for a 1.5× miss;
- r80 **1.50** at c = 1 %, down from 1.59, and 1.38 / 1.31 at c = 2 / 4 %.

The binding pool is decided in [Final pool: re-run the coverage proxy and apply the stop rule](https://github.com/dimadatascience/NanoGen/issues/86) on the regenerated census BAMs.

As in the census, the limit is biological. Few negatives (lymphocytes) are covered: 43 targets, 16 with ≥ 10 covered negatives, the best with 248.

## Method

**Tagged BAMs from FASTQ.** The driver is the one from [Regenerate census tagged BAMs with the short-read known barcode list](https://github.com/dimadatascience/NanoGen/issues/85): prototype `aa30310` (genotype-v3 `61d8809` + optional `known_bc`) plus its `regen.nf`. It runs BLAZE `--no-whitelisting --known-bc-list` → MERGE_READS → MINIMAP2 → FIXTAGS, then sort and index, and stops before genotyping. It uses the pinned image `sha256:764368d6…`.
- **Known list:** the sample's short-read QC list, as adopted in [Cell barcodes: short-read known list instead of BLAZE whitelisting](https://github.com/dimadatascience/NanoGen/issues/82). The restart samples are on the same footing as the regenerated census.
- **Lanes:** one FASTQ lane per sample.

**Mock first.** The first 250k reads per sample went through the whole driver and the coverage proxy:
- all 15 tasks exited 0;
- BLAZE assigned 17–26 % of reads to QC barcodes, so the QC lists belong to these FASTQs;
- peak memory was BLAZE ≤ 3.1 GB and MINIMAP2 ≤ 16.5 GB.

Full runs used the census memory settings (pilot peaks + 50 %) and times scaled from the mock rate.

| sample | reads | mean length | BLAZE | MINIMAP2 (peak RSS) | head wall time |
|---|---|---|---|---|---|
| AML4 | 53.1 M | 749 bp | 30 min | 13 min (21.4 GB) | 52 min |
| AML5 | 17.0 M | 811 bp | 10 min | 8 min (20.4 GB) | 24 min |
| sAML1 | 41.4 M | 902 bp | 27 min | 25 min (22.7 GB) | 63 min |

Every task completed on its first attempt.

**Power check.** `docs/research/power-check/power_check.py` at `9a2d94f` (byte-identical copy on the HPC), with the same parameters, verdict rule, c grid and seeds.
- Coverage is read from the BAM and no allele is read.
- Each sample is summarised raw, at worst case (0.685 / 0.709) and at central (0.9073 / 0.9513) scale.
- Tables come from `collect.py census`, unchanged.
- Pools (`fastq-restart/collect84.py`) sum the central-scale E_B over targets and take r80 from `power_check.power_numbers`. The census-only pool reproduces the power check exactly (E_B 22.4 / 37.1 / 52.4 and r80 1.59 / 1.45 / 1.39 at c = 1 / 2 / 4 %).

## Intake

| sample | QC cells | negatives (T/B/NK list, all in QC) | tumour (Leukemic) in QC | Ambiguous in QC | QC cells not annotated | targets |
|---|---|---|---|---|---|---|
| AML4 | 1,997 | 767 | **88** | 1,102 | 113 | 11 (11 SNV) |
| AML5 | 2,851 | 518 | 1,940 | 326 | 145 | 16 (14 SNV, 2 indel) |
| sAML1 | 4,292 | 1,073 | 2,438 | 281 | 810 | 16 (15 SNV, 1 indel) |

- **Sample names:** AML4, AML5 and sAML1, which equal the negatives CSV `sample` value.
- **Tumour compartment.** The census Seurat object doesn't hold these samples. Labels come from the earlier old-panel scRNA-seq study's TME object, where `compartment == "tumor"` is that study's tier-I malignant call. It is **phenotype-only**:
  - a myeloid, erythroid or HSC cell type;
  - in a cluster free of healthy-marrow cells;
  - expanded ≥ 4× over healthy marrow.

  No genotype enters the call, so it stays outcome-blind. `undefined` is kept apart, like Ambiguous.
- **Label conflicts.** 9 (AML5) and 28 (sAML1) tumour-labelled cells are in the provided negatives list, which comes from a newer annotation. The list wins. AML4's tumour compartment is thin: 88 cells, with most QC cells `undefined`. This matters for tumour-side criteria only.
- **Panel BEDs** have CRLF line endings. The proxy reads only the first five fields, but LF copies are used throughout.
- **Clonal status is unknown.** No bulk WES was found for these patients, so for now they are **calibration-only**, like PALM07E and PALM21E.

## Per-target power (central scale)

Verdicts: AML4 0/11 evaluable; AML5 0/16; sAML1 0/16, with 1 borderline. Targets with E_B > 0 are shown. The rest have ≤ 1 covered negative, or none.

| sample | target | covered negatives | B cells (median) | A UMIs (median) | E_B c = 0.25 % | 0.5 % | 1 % | 2 % | 4 % | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| AML4 | T1 | 110 | 51 | 108 | 0.25 | 0.49 | 0.94 | 1.48 | 2.04 | not evaluable |
| AML4 | T2 | 177 | 86 | 118 | 0.32 | 0.62 | 1.23 | 2.33 | 3.46 | not evaluable |
| AML4 | T3 | 5 | 3 | 4 | 0.01 | 0.02 | 0.04 | 0.07 | 0.11 | not evaluable |
| AML4 | T4 | 26 | 12 | 16 | 0.04 | 0.08 | 0.16 | 0.31 | 0.47 | not evaluable |
| AML4 | T6 | 9 | 4 | 5 | 0.01 | 0.03 | 0.05 | 0.11 | 0.18 | not evaluable |
| AML4 | T9 | 13 | 6 | 9 | 0.02 | 0.04 | 0.09 | 0.17 | 0.23 | not evaluable |
| AML4 | T10 | 1 | 1 | 0 | 0.00 | 0.01 | 0.02 | 0.04 | 0.07 | not evaluable |
| AML5 | T1 | 57 | 26 | 38 | 0.10 | 0.19 | 0.37 | 0.67 | 1.04 | not evaluable |
| AML5 | T2 | 2 | 1 | 5 | 0.00 | 0.01 | 0.02 | 0.04 | 0.07 | not evaluable |
| AML5 | T3 | 158 | 75 | 99 | 0.27 | 0.54 | 1.06 | 2.09 | 3.19 | not evaluable |
| AML5 | T4 | 11 | 5 | 14 | 0.04 | 0.07 | 0.15 | 0.22 | 0.14 | not evaluable |
| AML5 | T7 | 2 | 1 | 2 | 0.00 | 0.00 | 0.01 | 0.02 | 0.04 | not evaluable |
| AML5 | T10 | 1 | 0 | 0 | 0.00 | 0.00 | 0.00 | 0.01 | 0.02 | not evaluable |
| AML5 | T11 | 11 | 6 | 5 | 0.02 | 0.03 | 0.07 | 0.14 | 0.19 | not evaluable |
| AML5 | T16 | 58 | 28 | 60 | 0.16 | 0.31 | 0.60 | 0.76 | 1.01 | not evaluable |
| sAML1 | T1 | 42 | 21 | 28 | 0.07 | 0.15 | 0.29 | 0.57 | 0.80 | not evaluable |
| sAML1 | T2 | 6 | 3 | 3 | 0.01 | 0.01 | 0.03 | 0.05 | 0.11 | not evaluable |
| sAML1 | T4 | 2 | 1 | 1 | 0.00 | 0.00 | 0.01 | 0.02 | 0.04 | not evaluable |
| sAML1 | T7 | 47 | 22 | 34 | 0.10 | 0.20 | 0.30 | 0.55 | 0.95 | not evaluable |
| sAML1 | T8 | 82 | 38 | 68 | 0.17 | 0.34 | 0.61 | 1.16 | 1.59 | not evaluable |
| sAML1 | T9 | 248 | 118 | 185 | 0.47 | 0.92 | 1.82 | 3.38 | 4.84 | **borderline** |
| sAML1 | T10 | 54 | 27 | 42 | 0.11 | 0.21 | 0.41 | 0.67 | 1.17 | not evaluable |
| sAML1 | T11 | 107 | 50 | 104 | 0.24 | 0.43 | 0.84 | 1.43 | 2.02 | not evaluable |
| sAML1 | T12 | 28 | 13 | 30 | 0.06 | 0.11 | 0.19 | 0.36 | 0.49 | not evaluable |

"Covered negatives" is the raw proxy count; B cells, A UMIs and E_B are at central scale. The full columns (seeds, fail rates, c_min) are in `q84_restart/out/full/census.tsv`. sAML1 T9 passes the carried-over rule in 12/20 raw seeds and 0/20 at worst case, so it is borderline. Per-target evaluability is descriptive only (see [Underpowered census: stop, pool, or restart from FASTQ?](https://github.com/dimadatascience/NanoGen/issues/83)).

## Pools (E_B summed over targets, central scale)

| pool | E_B c = 0.25 % | 0.5 % | 1 % | 2 % | 4 % | r80 c = 0.25 % | 0.5 % | 1 % | 2 % | 4 % |
|---|---|---|---|---|---|---|---|---|---|---|
| AML4 | 0.66 | 1.28 | 2.51 | 4.51 | 6.56 | 6.53 | 4.30 | 3.15 | 2.53 | 2.26 |
| AML5 | 0.59 | 1.16 | 2.28 | 3.94 | 5.69 | 7.25 | 4.78 | 3.47 | 2.60 | 2.40 |
| sAML1 | 1.23 | 2.37 | 4.49 | 8.20 | 12.0 | 4.49 | 3.34 | 2.54 | 2.08 | 1.88 |
| restart (3) | 2.48 | 4.81 | 9.28 | 16.6 | 24.2 | 3.20 | 2.61 | 2.08 | 1.75 | 1.60 |
| census (5, canonical BAMs) | 6.36 | 12.2 | 22.4 | 37.1 | 52.4 | 2.33 | 1.86 | 1.59 | 1.45 | 1.39 |
| **census + restart (8), preview** | 8.83 | 17.0 | **31.7** | 53.7 | 76.7 | 2.06 | 1.72 | **1.50** | 1.38 | 1.31 |

The preview uses the census's canonical-BAM E_B. The regenerated census BAMs (×1.27 covered negatives on the PALM14E pilot) should raise it further, which "Final pool" will measure. The replacement stop rule (final-pool r80 > 2.0 at c = 1 %) is not close to triggering on this preview.

## HPC

Folder `q84_restart/` in the map's HPC working folder. Jobs:
- annotation 9064616; read counts 9064617; mock head 9064670 and coverage 9064736;
- full heads 9064762 (AML4), 9064763 (AML5), 9064764 (sAML1); coverage 9065622; collect 9065686.

Tagged BAMs: `runs/full_<S>/results/bams/<S>_knownbc_CB_UB.bam`. Labels: `derived/<S>_ct.csv`. LF BEDs: `derived/<S>.bed`.
