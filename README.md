
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




## 📤 Output Structure

The output structure:

```
params.outdir
|-- sample
|   |-- sample.csv
|   |-- logs
|       |-- blaze.log
|       |-- genotype.log
|       |-- target.log

```

- sample.csv contains the genotyping results for the targeted region. This CSV file includes the cell ID, the assigned genotype, the probability score of the genotype call, and the number of UMIs supporting the WT allele, the MUT allele, and unclassified reads.

- blaze.log: contains the statistics of the experiments, e.g. number of reads, high quality reads, number of cells and reads in polyT

- genotype.log: contains some statistics of genotyping, such as genotyped cells, median UMIs per cell

- target.log: contains the enrichment efficiency

## 🚀 Running the Pipeline

```bash
nextflow run path2/main.nf -c nextflow.config -profile singularity -entry LONG_READS  --input samples.csv --bedfile target.tsv
```

It requires the target file which is a tab-separated value file containing
chr, start position of target mutation, end position of target mutation, reference allele, alternative allele, gene associated

To be noted, start and end position are the same for SNP

## 👥 Authors

Andrea Cossa, Yinxiu Zhan

This pipeline was developed for high-throughput, single-cell transcriptomic analysis using Nanopore long reads and targeted region enrichment.