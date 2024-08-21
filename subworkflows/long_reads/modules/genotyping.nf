// GENOTYPING  module
 
nextflow.enable.dsl = 2
 
//


process GENOTYPING {
 
    tag "${sample}"
    publishDir "${params.outdir}/${sample}/", mode: 'copy'
 
    input:
    tuple val(sample), path(counts)
 
    output:
    tuple val(sample), path("${sample}.csv"), emit: genotypes
 
    script:
    """
    python ${baseDir}/bin/long_reads/genotyping.py ${counts} ${sample}
    """

    stub:
    """
    touch ${sample}.csv
    """
 
}