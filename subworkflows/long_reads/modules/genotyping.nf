// GENOTYPING  module
 
nextflow.enable.dsl = 2
 
//


process GENOTYPING {
 
    tag "${sample}"
 
    input:
    tuple val(sample), path(counts)
 
    output:
    tuple val(sample), path("${sample}.xlsx"), path("${sample}.tsv"), emit: genotypes
 
    script:
    """
    python ${baseDir}/bin/long_reads/genotyping.py ${counts} ${sample}
    """

    stub:
    """
    touch ${sample}.xlsx
    touch ${sample}.tsv
    """
 
}