// MINIMAP2 module

nextflow.enable.dsl = 2

//

process MINIMAP2 {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(fastq)

  output:
  tuple val(sample_name), path("unsorted.bam"), emit: bam

  script:
  """
  minimap2 \
  -ax splice \
  -t ${task.cpus} \
  ${params.minimap2_mindex} \
  ${fastq} \
  | samtools view -bS \
    -@ ${task.cpus} \
    -T ${params.minimap2_input_genome} \
    -o unsorted.bam
  """

  stub:
  """
  touch unsorted.bam
  """

}
