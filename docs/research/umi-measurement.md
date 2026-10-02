# UMI splitting and window double counting on real data (PALM14E)

Wayfinder ticket: *Measure UMI splitting and window double-counting on real data* (#9), map #3.
This is the measurement that the six checks in §4 of `umi-grouping-ont-10x.md` asked for (from the ticket *UMI grouping approaches for Nanopore 10x long reads*, #8).

**Patient data.** Everything below is aggregate. Targets appear as **T1–T7** in BED order (T3 is the 3-bp deletion; T3, T4 and T5 lie within 6 bp of each other). The mapping to coordinates, every per-cell table and every per-read table stays on the HPC in `PALM14E_subset/` (group `deo_p_dima_scm`).

## Setup

- **Input regenerated.** The original run's `work/` cache was gone, so `prep.nf` re-ran BLAZE → MERGE_READS → MINIMAP2 → FIXTAGS from `main`@`c442a32` with the PALM14E config. It ran in a separate launch, work and out dir (21.9 CPU-h).
- **Reusable subset.** Reads overlapping targets ±1 kb, with all flags kept: `PALM14E_subset/ontarget/palm14e_ontarget.bam`, 2.7 GB. It is reusable by the KIT-subset and regression tickets.
  - The full tagged BAM had 15.2 M primary records (matching the access note: 45 % of 33.9 M reads are in cells), 0 secondary and 0.76 M supplementary.
- **Stage 1 (`stage1_consensus.sh`).** Ran `main`'s own `make_consensus.sh` + `count_consensus.py` per cell with the baseline params: `nextflow.config` defaults `-q 20 -c 2 -w 100 -m 20 -d 1`, fmin 0.6. It kept the intermediates.
  - **Reproduction check:** across 4,947 cells, 97.1 % of cell × target rows are count-identical to the baseline `palm14e.csv`. UMI totals are WT 24,097 vs 24,051 and MUT 23,497 vs 23,497. Calls agree for 34,943 of 34,979 rows; the difference is BLAZE/umi_tools run-to-run noise.
- **Stage 2/3 (`measure_umi.py`, `followup.py`).** Re-group site-covering MAPQ ≥ 20 reads several ways (alleles from a pysam pileup that emulates bam-readcount), score each group with the uniform rule (n ≥ n_min, top allele ≥ fmin 0.6) and call genotypes with `main`'s `genotyping.py` logic, so only the grouping changes.
  - **Emulation check:** 98.7 % of group labels (92.7 % of exact counts) match bam-readcount over 65,408 groups.
  - The indel uses the tolerant ALT (k = 1) from #5.

Grouping schemes:

| name | key |
|---|---|
| current | BX × chr × ⌊POS/100⌋, all records (= `main`) |
| A | current, primary only (`-F 0x900`) |
| B1 | cell × target × exact raw UB, primary only |
| B1+C(d) | B1, then umi_tools-style directional clustering (n_a ≥ 2n_b − 1) within cell × target at Hamming ≤ 1 / LD ≤ 2 / LD ≤ 3 |

All UBs are 12 nt, so LD ≤ 1 is identical to Hamming ≤ 1.

## Findings by check

**1. umi_tools does almost nothing.** There are 214,141 bundles (cell × exact start position). 56 % hold one read and 70 % hold one UMI. Only **1.7 % of reads get BX ≠ UB**.

**2. Molecules straddle windows; one read end is fixed.**

- **Window splitting:** of exact-UB molecules with ≥ 2 reads, **31 % span > 1 window at 100 bp** (mean 1.44 windows) and 18 % do so at 1 kb.
- **By target:** the share ranges from 5 % (T2) to 60 % (T7). It is worse for reverse-strand reads (44 % vs 23 %).
- **Fixed read end:** within a molecule the rightmost aligned coordinate is fixed (median range 0–3 bp; 100 % of reads within ±5 bp). The leftmost POS, which the awk split keys on, varies (median range 6–15 bp, up to kb at T7).
- **UG:** a molecule with ≥ 2 reads carries a mean 6 distinct UG values, because umi_tools bundles by start.

**3. Raw UMIs form large error families.** Within cell × target, distinct-UB pairs at LD 1 occur at 0.09–0.95 % and at LD 2 at 0.8–4.7 %. In the cross-cell null they occur at 0.0005–0.002 % (LD 1) and 0.009–0.018 % (LD 2), which is **55–950× (LD 1) and 45–290× (LD 2) above the null**. Excess persists at LD 3 (30–90×).

The null rate is the false-merge rate. At LD ≤ 2 it is about 1–2 × 10⁻⁴ per molecule pair. With the clustered molecule load per cell × target (see check 5: median 1–3, p90 ≤ 9 molecules with ≥ 2 reads), that is below 1 % per cell × target.

**4. Merges are genuine siblings.**

- **Read-level purity:** reads absorbed from merged-in children disagree with their cluster's majority at **the same rate as reads inside one exact-UB molecule** (intrinsic). At LD ≤ 2:

  | target | T1 | T2 | T3 | T4 | T5 | T6 | T7 |
  |---|---|---|---|---|---|---|---|
  | intrinsic (%) | 4.9 | 0.06 | 8.8 | 7.7 | 6.5 | 12.7 | 11.1 |
  | merged children (%) | 5.4 | 0.08 | 9.9 | 9.0 | 7.7 | 12.4 | 11.4 |

  LD ≤ 3 adds slightly more (for example T6 13.1 %).
- **Parent–child concordance:** label discordance for children with ≥ 2 reads is 0–8 %, flat across child–parent distance 1–4. That is the level the intrinsic mixing predicts.
- **Earlier pairwise check:** sibling pairs (LD 1–3) were discordant at about half the random-pair rate, but that used n = 1 labels. It is explained by single-read error, not false merges.

**4b. Within-molecule allele mixing (new).** At het SNV targets, 5–13 % of reads inside one exact-UB molecule carry the other allele, against 0.06 % at T2 (near-clonal). That is far above ONT base error, so it points to PCR chimeras / template switching or UMI collisions. This matters for fmin. It also means **small fragments of a WT molecule can come out MUT-majority**, which is how `main` produces spurious single-UMI MUT calls.

**5. What n_min removes; effect on counts and calls.**

- **Counts are inflated.** `main` counts each LD ≤ 2 molecule **5.4–8.1× on average** (median 1, heavy tail), and **18–51 % of molecules are counted more than once**, through window splits plus error-UMI children that reach 2 reads.
- **Molecule totals:** cell × target molecules with ≥ 2 reads fall from about 62 k (B1) to about 13 k (B1+C LD ≤ 2).
- **Molecules per cell × target:** before clustering the maximum is 2,010 (median 1–14). After LD ≤ 2 clustering the median is 1–3, the p90 is 2–9 and the maximum is 47.
- **Reads per molecule:** median 3–4 reads; p90 up to 1,644.
- **Singletons after clustering** are 77–96 % of groups but hold **only 1.0–1.4 % of reads**. So n_min mostly removes unclusterable stray reads; n_min 2 vs 3 changes the clustered UMI totals by about 30 %.

UMI totals summed over targets (WT / MUT):

| scheme | n_min 1 | n_min 2 | n_min 3 |
|---|---|---|---|
| current | 113,360 / 110,629 | 30,534 / 29,898 | 19,391 / 18,977 |
| B1 | 109,314 / 106,690 | 28,674 / 28,618 | 18,223 / 18,105 |
| B1+C Hamming ≤ 1 | 62,802 / 59,854 | 17,045 / 16,326 | 11,925 / 11,565 |
| B1+C LD ≤ 2 | 32,090 / 30,683 | 6,145 / 5,550 | 4,240 / 3,808 |
| B1+C LD ≤ 3 | 22,988 / 22,566 | 4,075 / 3,666 | 3,079 / 2,684 |

Calls, current (n_min 2, = baseline) → B1+C LD ≤ 2 (n_min 2), over 34.9 k cell × target rows:

| current → new | MUT | NA | WT |
|---|---|---|---|
| **MUT** (2,819) | 1,863 | 652 | 304 |
| **WT** (647) | 37 | 235 | 375 |
| **NA** (31,016) | 171 | 30,825 | 20 |

- **MUT → NA/WT:** about a third of current MUT calls change, consistent with spurious MUT "UMIs" from fragmented WT molecules (4b) and inflated counts. The median WT+MUT molecule count of a covered cell × target drops to 0.29× of current.
- **T2 (near-clonal) is stable**, while het targets move most.
- **Caveat:** `main`'s caller calls MUT on a single MUT UMI, so it amplifies any grouping noise. The spec's per-target error and binomial tail will behave differently. Re-score with the new caller in the regression check.

**6. Supplementary records act as reads.** There are no secondary records (minimap2 `-ax splice` output); supplementary records are 1.1 % of MAPQ ≥ 20 records (3.3 % at T7). A read plus its supplementary share UB and window, so **a single read passes n_min = 2**. Option A alone (`-F 0x900`) halves T2's MUT calls (562 → 284) and moves 506 current MUT calls to NA. The filter is needed whatever the grouping.

## What this decides for the grouping redesign

- **Option A (filter `-F 0x900`, drop `--paired`)** is required: supplementary records inflate n_min.
- **Option B (cell × target key instead of start/window)** is supported by checks 1 and 2: umi_tools bundles are 56 % singletons, and 31 % of molecules straddle 100 bp windows. If a position term is kept, the fixed right end is the stable one.
- **Option C (directional clustering)** is strongly supported (checks 3–4). **LD ≤ 2** captures most of the excess with no purity loss. LD ≤ 3 adds a little more merging at slightly lower purity, and Hamming ≤ 1 leaves about 2× more fragments. The estimated false-merge risk is below 1 % per cell × target.
- **n_min after clustering** removes about 1 % of reads. A sub-n_min "low-confidence" class would mostly hold stray reads.
- **fmin needs a look** (4b): within-molecule mixing of 5–13 % at het targets.

## Reproduce (HPC, group `deo_p_dima_scm`)

```bash
cd /hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E_subset
sbatch run_prep.sh                    # NanoGen/prep.nf at c442a32 -> ontarget/
sbatch scripts/stage1_consensus.sh    # per-cell main consensus, intermediates kept in cells/
sbatch scripts/stage2_measure.sh      # checks 1-6 -> measure_out/
sbatch scripts/stage3_followup.sh     # F1-F4 -> measure_out/ (caches records.pkl)
```

Scripts are in `docs/research/umi-measurement/`. `synth_test.py` builds a synthetic fixture for smoke-testing.
