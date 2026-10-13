
# 📡 Long Reads Single-Cell RNA-seq Pipeline (Nanopore)

This **Nextflow pipeline** performs an end-to-end analysis of **single-cell RNA sequencing (scRNA-seq) data generated with Oxford Nanopore Technologies**, focused on targeted/enriched transcript regions. The workflow includes demultiplexing, alignment, consensus generation, and genotyping.

## 🔄 Workflow Overview

```nextflow
workflow long_reads {
    take:
        ch_input 

    main:
        BLAZE → MERGE_READS → MINIMAP2 → SPLIT_BAM 
              → CONSENSUS → GENOTYPING

    emit:
        results = GENOTYPING.out.genotypes
}
```

## 🧬 Pipeline Steps

1. **BLAZE**: Barcode-aware filtering and demultiplexing of long reads.
2. **MERGE_READS**: Merges filtered reads per sample.
3. **MINIMAP2**: Aligns reads to the reference genome using `minimap2`.
5. **SPLIT_BAM**: Splits BAMs by cell barcode.
6. **CONSENSUS**: Builds consensus transcripts for each single cell.
7. **GENOTYPING**: Derives genotypes from allele count data.

## 📥 Input

The pipeline expects an input file with header

| sample         | lane              | fastq          |
| -------------- | ----------------- | -------------- |
| sample1        | l1                | path2/fastq.gz |
| .....          | .....             | .....          |

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
| `P_s` | Error rate of the target used for the call, as in `sample_targets.csv` |
| `p` | Binomial upper tail P(X ≥ MUT), X ~ Bin(MUT + WT, `P_s`): the probability of at least this many MUT UMIs if the cell is WT. Empty (NaN) when MUT + WT = 0 |
| `genotype` | Call: `MUT` if MUT ≥ `min_mut_umis` and p < 0.1; `WT` if MUT = 0 and WT ≥ 2, or MUT ≥ `min_mut_umis` and p ≥ 0.9; empty (NA) when undetermined |
| `het_detect_lb` | For WT calls, 1 − 0.5^WT: a lower bound on the probability that a heterozygous cell would have shown a MUT UMI. Empty otherwise |

- sample_targets.csv contains the error estimate of each target, pooled over cells, one row per target, written without an index column:

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

- blaze.log: contains the statistics of the experiments, e.g. number of reads, high quality reads, number of cells and reads in polyT

- genotype.log: contains some statistics of genotyping, such as genotyped cells, median UMIs per cell

- target.log: contains the enrichment efficiency

## ⚙️ Parameters

`nextflow_schema.json` is the source of truth for parameters. They are validated before any process runs, so a mistyped value fails in seconds. The visible parameters are:

| Parameter | Default | Spec symbol | Meaning |
| --------- | ------- | ----------- | ------- |
| `--input` | — | — | Samplesheet CSV: sample, lane, fastq; or sample, bam with `--input_type tagged_bam` |
| `--input_type` | `fastq` | — | `fastq`, or `tagged_bam` to start from tagged, aligned BAMs at SPLIT_BAM (see Input) |
| `--outdir` | `results` | — | Directory where results are saved |
| `--minimap2_input_genome` | — | — | Reference genome FASTA |
| `--minimap2_mindex` | — | — | minimap2 index of the reference genome; required with `--input_type fastq` |
| `--bedfile` | — | — | Tab-separated target file: chr, start, end, ref, alt, gene |
| `--min_mapping_quality` | `20` | — | Minimum MAPQ for a primary mapped read to enter UMI grouping |
| `--min_reads` | `3` | n_min | Minimum usable reads supporting an allele to label a UMI cluster |
| `--min_fraction` | `0.6` | f_min | Minimum fraction of usable reads supporting an allele to label a UMI cluster, in (0.5, 1] |
| `--snv_min_base_quality` | `20` | — | Minimum base quality (Phred) for a read to be usable at an SNV target |
| `--error_floor` | `0.0001` | ε | Lower bound on the per-target error rate, in (0, 1) |
| `--min_mut_umis` | `1` | m_min | Minimum number of MUT UMIs for a MUT call; set `2` to check results don't rest on single-UMI MUT calls |

Hidden parameters are listed by `--help --show_hidden_params`. They include the UMI grouping rule's thresholds `--umi_max_ld` (`2`), `--umi_ratio` (`2`) and `--umi_small_reads` (`3`), and the BLAZE options `--blaze_expect_cells` and `--blaze_edit_distance`.

Retired parameters still run but log a warning and are ignored: `--baseline_error` is now `--error_floor`, `--min_base_quality` is now `--snv_min_base_quality`, and `--window_width` and `--min_dist_umi_collapsing` are no longer parameters.

## 🚀 Running the Pipeline

```bash
nextflow run path2/main.nf -c nextflow.config -profile singularity -entry LONG_READS  --input samples.csv --bedfile target.tsv
```

It requires the target file which is a tab-separated value file containing
chr, start position of target mutation, end position of target mutation, reference allele, alternative allele, gene associated

To be noted, start and end position are the same for SNP

### Container

Every process runs in one pinned image, `ghcr.io/dimadatascience/nanogen`, referenced by digest through `params.container` in `nextflow.config`. Use `-profile docker` or `-profile singularity`. The image is built from `docker/` by the `docker` GitHub workflow, which reports the new digest; bump `params.container` (and its default in `nextflow_schema.json`) in a PR to adopt it.

## 👥 Authors

Andrea Cossa, Yinxiu Zhan

This pipeline was developed for high-throughput, single-cell transcriptomic analysis using Nanopore long reads and targeted region enrichment.