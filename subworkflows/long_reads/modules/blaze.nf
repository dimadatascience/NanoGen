// BLAZE module

nextflow.enable.dsl = 2

//

process BLAZE {

  tag "${sample_name}: ${lane}"

  input:
  tuple val(sample_name), val(lane), path(fastq)

  output:
  tuple val(sample_name), path("${sample_name}_${lane}_matched_reads.fastq.gz"), emit: filtered_fastq

  script:
  """
  # Default option... --expect-cells #params.blaze_expect_cells
  blaze ${fastq} \
  --no-whitelisting \
  --threads ${task.cpus} \
  --output-prefix ${sample_name}_${lane}_ \
  --max-edit-distance ${params.blaze_edit_distance}
  """

  stub:
  """
  touch ${sample_name}_${lane}_matched_reads.fastq.gz
  """

}