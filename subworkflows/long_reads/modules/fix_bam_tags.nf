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
  jvarkit samjdk \
  -e 'String s=record.getReadName(); int h=s.indexOf("#"),u=s.indexOf("_");record.setReadName(s.substring(h+1));record.setAttribute("CB",s.substring(0,u));record.setAttribute("UB",s.substring(u+1,h));return record;' \
  --samoutputformat BAM \
  -o CB_UB.bam \
  ${bam}
  """

  stub:
  """
  touch CB_UB.bam
  """

}
