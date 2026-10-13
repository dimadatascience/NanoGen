// GENOTYPING  module
 
nextflow.enable.dsl = 2
 
//


process GENOTYPING {
 
    tag "${sample}"
    publishDir "${params.outdir}/${sample}/", mode: 'copy', pattern: "${sample}*.csv"
    publishDir "${params.outdir}/${sample}/logs/", mode: 'copy', pattern: "genotype.log"
 
    input:
    tuple val(sample), path(counts)
 
    output:
    tuple val(sample), path("${sample}.csv"), emit: genotypes
    tuple val(sample), path("${sample}_targets.csv"), emit: targets
    tuple val(sample), path("genotype.log"), emit: log
 
    script:
    """
    python ${baseDir}/bin/long_reads/genotyping.py --input ${counts} --output ${sample}.csv --targets_output ${sample}_targets.csv --error_floor ${params.error_floor} --min_mut_umis ${params.min_mut_umis} --snv_min_base_quality ${params.snv_min_base_quality}
    """

    stub:
    """
    touch ${sample}.csv
    touch ${sample}_targets.csv
    touch genotype.log
    """
 
}