// MERGE_READS module

nextflow.enable.dsl = 2

//

process MERGE_READS {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(fastqs)

  output:
  tuple val(sample_name), path("matched_reads.fastq.gz"), emit: filtered_fastq

  script:
  """
  zcat ${fastqs} | pigz --fast > matched_reads.fastq.gz
  """

  stub:
  """
  touch matched_reads.fastq.gz
  """

}
