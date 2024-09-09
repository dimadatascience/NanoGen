// BLAZE module

nextflow.enable.dsl = 2

//

process BLAZE {

  tag "${sample_name}: ${lane}"
  publishDir "${params.outdir}/${sample_name}/logs/", mode: 'copy', pattern: 'blaze.log'

  input:
  tuple val(sample_name), val(lane), path(fastq)

  output:
  tuple val(sample_name), path("${sample_name}_${lane}_matched_reads.fastq.gz"), emit: filtered_fastq
  tuple val(sample_name), path("blaze.log"), emit: log

  script:
  """
  blaze \
  --expect-cells ${params.blaze_expect_cells} \
  --threads ${task.cpus} \
  --output-prefix ${sample_name}_${lane}_ \
  --max-edit-distance ${params.blaze_edit_distance} \
  ${fastq}

  cp .command.out blaze.log
  """

  stub:
  """
  touch ${sample_name}_${lane}_matched_reads.fastq.gz
  """

}
