// FIX_BAM_TAGS module

nextflow.enable.dsl = 2

//

process FIXTAGS {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(bam)

  output:
  tuple val(sample_name), path("CB_UB.bam"), emit: bam

  script:
  """
  python ${baseDir}/bin/long_reads/fix_tags.py ${bam} 
  """

  stub:
  """
  touch CB_UB.bam
  """

}
