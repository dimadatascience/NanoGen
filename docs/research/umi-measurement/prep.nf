// Regenerate the CB/UB-tagged BAM for PALM14E with main@c442a32 modules,
// then keep only reads overlapping targets +-1 kb (all flags kept).
nextflow.enable.dsl = 2
include { BLAZE } from "./subworkflows/long_reads/modules/blaze.nf"
include { MERGE_READS } from "./subworkflows/long_reads/modules/merge_reads.nf"
include { MINIMAP2 } from "./subworkflows/long_reads/modules/minimap2.nf"
include { FIXTAGS } from "./subworkflows/long_reads/modules/fix_bam_tags.nf"

process ONTARGET {
  publishDir "${params.outdir}", mode: 'copy'
  cpus 8
  memory '20 GB'
  time '3h'
  container 'docker://acox1/scmseq:1.0'
  input:
  tuple val(sample_name), path(bam)
  output:
  path("${sample_name}_ontarget.bam*")
  path("${sample_name}_flagstat.txt")
  script:
  """
  awk 'BEGIN{OFS="\\t"}{s=\$2-1000; if(s<0)s=0; print \$1,s,\$3+1000}' ${params.bedfile} > padded.bed
  samtools sort -@ ${task.cpus} -m 2G -o sorted.bam ${bam}
  samtools index sorted.bam
  samtools flagstat sorted.bam > ${sample_name}_flagstat.txt
  samtools view -@ ${task.cpus} -b -L padded.bed -o ${sample_name}_ontarget.bam sorted.bam
  samtools index ${sample_name}_ontarget.bam
  rm sorted.bam sorted.bam.bai
  """
}

workflow {
  ch = Channel.fromPath(params.input).splitCsv(header: true)
        .map { row -> tuple(row.sample, row.lane, file(row.fastq)) }
  BLAZE(ch)
  MERGE_READS(BLAZE.out.filtered_fastq.groupTuple(by: 0))
  MINIMAP2(MERGE_READS.out.filtered_fastq)
  FIXTAGS(MINIMAP2.out.bam)
  ONTARGET(FIXTAGS.out.bam)
}
