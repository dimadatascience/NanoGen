// GENOTYPING  module
 
nextflow.enable.dsl = 2
 
//


process GENOTYPING {
 
    tag "${sample}"
    publishDir "${params.outdir}/${sample}/", mode: 'copy', pattern: "${sample}.csv"
    publishDir "${params.outdir}/${sample}/logs/", mode: 'copy', pattern: "genotype.log"
 
    input:
    tuple val(sample), path(counts)
 
    output:
    tuple val(sample), path("${sample}.csv"), emit: genotypes
    tuple val(sample), path("genotype.log"), emit: log
 
    script:
    """
    python ${baseDir}/bin/long_reads/genotyping.py --input ${counts} --output ${sample}.csv --baseline_error ${params.baseline_error}
    """

    stub:
    """
    touch ${sample}.csv
    touch genotype.log
    """
 
}