
# 📡 Long Reads Single-Cell RNA-seq Pipeline (Nanopore)

This **Nextflow pipeline** performs an end-to-end analysis of **single-cell RNA sequencing (scRNA-seq) data generated with Oxford Nanopore Technologies**, focused on targeted/enriched transcript regions. The workflow includes demultiplexing, alignment, consensus generation, and genotyping.

## 🔄 Workflow Overview

```nextflow
workflow long_reads {
    take:
        ch_input 

    main:
        BLAZE → MERGE_READS → MINIMAP2 → FIXTAGS → SPLIT_BAM
              → CONSENSUS → COLLAPSE → GENOTYPING

    emit:
        results = GENOTYPING.out.genotypes
}
```

## 🧬 Pipeline Steps

1. **BLAZE**: Barcode-aware filtering and demultiplexing of long reads.
2. **MERGE_READS**: Merges filtered reads per sample.
3. **MINIMAP2**: Aligns reads to the reference genome using `minimap2`.
4. **FIXTAGS**: Writes the cell barcode and UMI into the `CB` and `UB` tags.
5. **SPLIT_BAM**: Splits BAMs by cell barcode and reports on-target reads.
6. **CONSENSUS**: Groups reads per cell × target into consensus UMIs and labels each WT, MUT or OTH.
7. **COLLAPSE**: Gathers the consensus UMI counts of all cells of a sample.
8. **GENOTYPING**: Estimates the error rate and the contamination rate of each target and calls each cell × target MUT, WT or undetermined.

## 🧮 Genotyping model

The model is specified in [`docs/genotyping/genotyping.pdf`](docs/genotyping/genotyping.pdf) (source `genotyping.tex`); this is a summary.
Within each cell, the reads covering a target are grouped by UMI, merging UMIs that differ only by sequencing errors, so each group stands for one molecule.
A group with enough reads, most of them agreeing, becomes a consensus UMI labelled WT (reference allele), MUT (alternative allele) or OTH (any other allele); groups with too few reads are counted but not used.
OTH UMIs can only arise from errors, so the error rate `P_s` of each target is estimated from its OTH UMIs pooled over all cells, with a lower floor.
A MUT UMI in a cell without the mutation can also come from contamination (ambient RNA, a PCR chimera, barcode swapping), which gives a clean MUT consensus that `P_s` cannot see. With `--negative_cells`, a list of cells that cannot carry the mutations (T and B cells), the contamination rate `c_s` of each target is the MUT fraction of those cells' WT+MUT UMIs; without the file `c_s` = 0.
For each cell × target, a binomial tail with `P_null` = `P_s` + `c_s` gives the probability of seeing at least the observed MUT UMIs if the cell were WT.
A cell is called MUT when that probability is small, WT when it has no MUT UMI and at least two WT UMIs or when the MUT UMIs are well explained by error and contamination, and is left undetermined otherwise.
So a single MUT UMI is accepted in a cell with few WT UMIs and rejected in a cell with many, without a cell-type-specific threshold: at `c_s` = 0.03, one MUT UMI next to 0–2 WT UMIs is MUT and next to 3 or more is undetermined.
Single-UMI MUT calls, MUT calls left in negative-control cells and the chance of missing a heterozygous cell are reported alongside the calls.

## 📥 Input

The pipeline expects an input file with header

| sample         | lane              | fastq          |
| -------------- | ----------------- | -------------- |
| sample1        | l1                | path2/fastq.gz |
| .....          | .....             | .....          |

An optional `known_bc` column gives the sample's short-read cell barcodes (for example Cell Ranger's filtered `barcodes.tsv.gz`, or a QC-passed list; plain or gzipped, one barcode per line, the `-1` suffix is stripped). With it, BLAZE assigns reads to that list (`--no-whitelisting --known-bc-list`) instead of building its own whitelist from the long-read barcode counts, so cells with few reads at the targets are not dropped by the count cutoff. Give the same file on every lane of a sample, or leave it empty.

### Starting from a tagged BAM

With `--input_type tagged_bam` the run skips BLAZE, MERGE_READS, MINIMAP2 and FIXTAGS and starts at SPLIT_BAM. The samplesheet has one row per sample:

| sample         | bam                  |
| -------------- | -------------------- |
| sample1        | /path2/sample1.bam   |
| .....          | .....                |

Each BAM holds one sample's reads aligned to `--minimap2_input_genome`, with the cell barcode in the `CB` tag and the UMI in the `UB` tag, as FIXTAGS emits it (`CB_UB.bam` in its work directory). It need not be sorted or indexed. `--minimap2_mindex` is not needed, and `blaze.log` is not written. Two runs from the same tagged BAM see identical reads, without BLAZE run-to-run noise.



## 📤 Output Structure

The output structure:

```
params.outdir
|-- sample
|   |-- sample.csv
|   |-- sample_targets.csv
|   |-- logs
|       |-- blaze.log
|       |-- genotype.log
|       |-- target.log
```

- sample.csv contains the genotyping results for the targeted region, one row per cell × target, written without an index column:

| Column | Meaning |
| ------ | ------- |
| `cell` | Cell barcode |
| `target_id` | Target as `chr:start:ref>alt`: the 1-based bam-readcount position, with `ref` and `alt` as written in the BED (e.g. `chr4:54733160:GAT>-GAT`) |
| `gene` | Gene name from the BED |
| `type` | `SNV`, or `INDEL` for a `+`/`-` alt |
| `WT` | Number of UMIs with consensus label WT |
| `MUT` | Number of UMIs with consensus label MUT |
| `OTH` | Number of UMIs with consensus label OTH (third allele) |
| `P_s` | Error rate of the target, as in `sample_targets.csv` |
| `c_s` | Contamination rate of the target, as in `sample_targets.csv`; 0 without `--negative_cells` |
| `P_null` | `P_s` + `c_s`, capped at 1: the probability that a UMI of a cell without the mutation reads MUT, used for the call |
| `p` | Binomial upper tail P(X ≥ MUT), X ~ Bin(MUT + WT, `P_null`): the probability of at least this many MUT UMIs if the cell is WT. Empty (NaN) when MUT + WT = 0 |
| `genotype` | Call: `MUT` if MUT ≥ `min_mut_umis` and p < 0.1; `WT` if MUT = 0 and WT ≥ 2, or MUT ≥ `min_mut_umis` and p ≥ 0.9; empty (NA) when undetermined |
| `het_detect_lb` | For WT calls, 1 − 0.5^WT: a lower bound on the probability that a heterozygous cell would have shown a MUT UMI. Empty otherwise |
| `negative_control` | `True` for a cell listed in `--negative_cells`; such cells are called like any other |

- sample_targets.csv contains the error estimate and the limitations of each target, pooled over cells, one row per target, written without an index column:

| Column | Meaning |
| ------ | ------- |
| `target_id` | Target, as in `sample.csv` |
| `gene` | Gene name from the BED |
| `type` | `SNV` or `INDEL` |
| `N` | Labelled (WT, MUT or OTH) UMIs over all cells; discarded and sub-n_min clusters are not counted |
| `w` | OTH UMIs over all cells |
| `d` | Error channels: 2 for an SNV, 1 for an indel |
| `q_hat` | Estimated error rate w / (d · N); empty (NaN) when N = 0 |
| `eps` | Error floor, `--error_floor` |
| `P_s` | Error rate used for the calls, max(`q_hat`, `eps`) |
| `floor_applied` | `True` when `P_s` is the floor, i.e. `q_hat` < `eps` or N = 0 |
| `n_neg_cells` | Negative-control cells with at least one WT or MUT UMI at the target |
| `n_neg_umis` | WT + MUT UMIs of the negative-control cells at the target |
| `n_neg_mut` | MUT UMIs of the negative-control cells at the target |
| `c_s` | Contamination rate `n_neg_mut` / `n_neg_umis`; 0 when `n_neg_umis` = 0, including without `--negative_cells`. A 0 from few negative UMIs is a weak estimate, not evidence of no contamination |
| `P_null` | `P_s` + `c_s`, capped at 1 |
| `min_base_qual` | Base-quality filter applied to the target's reads: `--snv_min_base_quality` at an SNV, 0 at an indel |
| `wt_reads` | Usable reads of the WT-labelled clusters over all cells. Indel diagnostic: empty (NaN) for an SNV |
| `wt_altlike_reads` | ALT-like reads: of `wt_reads`, the reads matching the alt by the tolerant match, showing WT→ALT leakage at a repeat. Empty for an SNV |
| `wt_altlike_frac` | `wt_altlike_reads` / `wt_reads`; empty for an SNV or when `wt_reads` = 0 |
| `mixed_discarded` | Mixed discarded clusters: discarded clusters (≥ n_min usable reads, no allele reaching n_min and f_min) with both REF and ALT usable reads |
| `subnmin_umis` | Sub-n_min clusters: clusters covering the target with fewer than n_min usable reads, never labelled |
| `subnmin_reads` | Usable reads of the sub-n_min clusters |
| `n_MUT` | Cells called MUT |
| `n_WT` | Cells called WT |
| `n_undet` | Cells left undetermined, including cells with no UMI at the target |
| `n_MUT_single_umi` | MUT calls resting on a single MUT UMI (MUT = 1), the calls most exposed to contamination |
| `n_MUT_neg` | MUT calls in negative-control cells: contamination the caller still lets through |
| `one_umi_mut_frac_neg` | The one-UMI test: among negative-control cells with exactly one WT or MUT UMI at the target, the fraction whose UMI is MUT (an independent view of `c_s`). Empty without such cells |
| `one_umi_mut_frac_other` | The same fraction among the other cells. A single MUT UMI carries information at this target only where this is well above `one_umi_mut_frac_neg`; when the two are close, single-UMI calls are contamination whatever the threshold |

- blaze.log: contains the statistics of the experiments, e.g. number of reads, high quality reads, number of cells and reads in polyT

- genotype.log: contains statistics of genotyping. The first lines report the number of cells, the total and median UMIs per cell, the median number of genes covered per cell and the fraction of cells with UMIs. The next lines summarise the calls and their limitations:
  - the MUT, WT and undetermined call totals
  - the number of MUT calls with a single MUT UMI (m = 1), and their fraction of MUT calls
  - the number of WT calls with k = 2, k = 3 and k ≥ 4 WT UMIs
  - the median and minimum `het_detect_lb` over WT calls
  - the targets at the error floor ε, their number out of all targets and their IDs
  - the `min_mut_umis` value used
  - the number of negative-control cells, the targets with `c_s` > 0 and their `c_s`, and the MUT calls left in negative-control cells

  A statistic over no calls is reported as `nan`.

- target.log: contains the enrichment efficiency

## ⚙️ Parameters

`nextflow_schema.json` is the source of truth for parameters. They are validated before any process runs, so a mistyped value fails in seconds. The visible parameters are:

| Parameter | Default | Spec symbol | Meaning |
| --------- | ------- | ----------- | ------- |
| `--input` | — | — | Samplesheet CSV: sample, lane, fastq (optional known_bc); or sample, bam with `--input_type tagged_bam` |
| `--input_type` | `fastq` | — | `fastq`, or `tagged_bam` to start from tagged, aligned BAMs at SPLIT_BAM (see Input) |
| `--outdir` | `results` | — | Directory where results are saved |
| `--minimap2_input_genome` | — | — | Reference genome FASTA |
| `--minimap2_mindex` | — | — | minimap2 index of the reference genome; required with `--input_type fastq` |
| `--bedfile` | — | — | Tab-separated target file: chr, start, end, ref, alt, gene |
| `--negative_cells` | — | c_s | CSV of negative-control cells (sample, cell), e.g. T and B cells; their MUT UMIs give the contamination rate c_s per target. Omit it to keep c_s = 0 |
| `--min_mapping_quality` | `20` | — | Minimum MAPQ for a primary mapped read to enter UMI grouping |
| `--min_reads` | `3` | n_min | Minimum usable reads supporting an allele to label a UMI cluster |
| `--min_fraction` | `0.6` | f_min | Minimum fraction of usable reads supporting an allele to label a UMI cluster, in (0.5, 1] |
| `--snv_min_base_quality` | `20` | — | Minimum base quality (Phred) for a read to be usable at an SNV target |
| `--error_floor` | `0.0001` | ε | Lower bound on the per-target error rate, in (0, 1) |
| `--min_mut_umis` | `1` | m_min | Minimum number of MUT UMIs for a MUT call; set `2` to check results don't rest on single-UMI MUT calls |

Hidden parameters are listed by `--help --show_hidden_params`. They include the UMI grouping rule's thresholds `--umi_max_ld` (`2`), `--umi_ratio` (`2`) and `--umi_small_reads` (`3`), and the BLAZE options `--blaze_expect_cells` and `--blaze_edit_distance`.

Retired parameters still run but log a warning and are ignored: `--baseline_error` is now `--error_floor`, `--min_base_quality` is now `--snv_min_base_quality`, and `--window_width` and `--min_dist_umi_collapsing` are no longer parameters.

## 🔁 Changes in genotype-v3

genotype-v3 adds contamination to the caller's null and changes nothing else:

- **`--negative_cells`.** A CSV (`sample`, `cell`) of cells that cannot carry the targeted mutations, typically T and B cells from the short-read annotation. Per target, the MUT fraction of their WT+MUT UMIs is the contamination rate `c_s`, and the call uses `P_null` = `P_s` + `c_s` instead of `P_s`. Without the file `c_s` = 0 and the calls are identical to genotype-v2.
- **Why.** `P_s` is estimated from third-allele UMIs and sits near its floor; a contaminating molecule gives a clean MUT consensus, so a single MUT UMI was always significant, in T cells as in blasts. Requiring two MUT UMIs removed the T-cell calls but also most blast calls, because at 2–5 UMIs per cell a heterozygous cell shows two MUT UMIs less than half the time. With `c_s` in the null, a lone MUT UMI is accepted when the cell has few WT UMIs and rejected when it has many.
- **New columns.** `c_s`, `P_null` and `negative_control` in `sample.csv`; `n_neg_cells`, `n_neg_umis`, `n_neg_mut`, `c_s`, `P_null`, `n_MUT_neg`, `one_umi_mut_frac_neg` and `one_umi_mut_frac_other` in `sample_targets.csv`; three lines at the end of `genotype.log`. Column positions after `P_s` shift accordingly.
- **The one-UMI test.** `one_umi_mut_frac_other` against `one_umi_mut_frac_neg` says, per patient and target, whether single-UMI calls carry information before any threshold is chosen.

## 🔁 Changes in genotype-v2

For users of `main`, genotype-v2 changes the outputs, the parameters and the calls:

- **MIS → OTH.** The third-allele UMI count column is renamed from `MIS` to `OTH`.
- **`target_id`.** A target is identified by one `target_id` column (`chr:start:ref>alt`) instead of `chr`, `start`, `end`, `ref` and `alt`, next to `gene` and a new `type`; `sample.csv` no longer has an index column.
- **New `sample_targets.csv`.** The error rate is estimated per target, pooled over cells, instead of per cell, pooled over targets; this file reports it with each target's limitations (see Output Structure).
- **Parameters.** `--baseline_error` is now `--error_floor`, `--min_base_quality` is now `--snv_min_base_quality`, and `--window_width` and `--min_dist_umi_collapsing` are retired; `--min_mut_umis` is new. Retired names still run with a warning. See Parameters.
- **Tagged-BAM entry point.** `--input_type tagged_bam` starts a run from aligned BAMs with `CB`/`UB` tags (see Input).
- **Calls differ from `main`.** UMIs are grouped per cell × target instead of per 100 bp window, so a molecule is counted once; indel alleles use the same consensus rule as SNVs; the error rate is per target. On PALM14E (4,947 cells, both runs from the same tagged BAM), MUT calls drop from 2,289 to 1,393 and 1,073 cells change at least one call; 174 calls flip MUT → WT and 5 WT → MUT. Each change is attributed to its cause, and genotype-v2 is checked against the spec, in [`docs/research/regression-palm14e.md`](docs/research/regression-palm14e.md).

## 🚀 Running the Pipeline

```bash
nextflow run path2/main.nf -c nextflow.config -profile singularity -entry LONG_READS  --input samples.csv --bedfile target.tsv
```

To model contamination, add `--negative_cells negative_cells.csv`, a CSV with a header `sample,cell` listing, per sample, the barcodes of cells that cannot carry the targeted mutations (T and B cells). Samples absent from the file run with `c_s` = 0.

It requires the target file which is a tab-separated value file containing
chr, start position of target mutation, end position of target mutation, reference allele, alternative allele, gene associated

To be noted, start and end position are the same for SNP

### Container

Every process runs in one pinned image, `ghcr.io/dimadatascience/nanogen`, referenced by digest through `params.container` in `nextflow.config`. Use `-profile docker` or `-profile singularity`. The image is built from `docker/` by the `docker` GitHub workflow, which reports the new digest; bump `params.container` (and its default in `nextflow_schema.json`) in a PR to adopt it.

## 👥 Authors

Andrea Cossa, Yinxiu Zhan

This pipeline was developed for high-throughput, single-cell transcriptomic analysis using Nanopore long reads and targeted region enrichment.