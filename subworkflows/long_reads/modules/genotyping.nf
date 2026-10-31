// GENOTYPING  module
 
nextflow.enable.dsl = 2
 
//


process GENOTYPING {
 
    tag "${sample}"
    publishDir "${params.outdir}/${sample}/", mode: 'copy', pattern: "${sample}*.csv"
    publishDir "${params.outdir}/${sample}/logs/", mode: 'copy', pattern: "genotype.log"
 
    input:
    tuple val(sample), path(counts)
    path negative_cells  // CSV (sample, cell) of negative-control cells, or [] for none
 
    output:
    tuple val(sample), path("${sample}.csv"), emit: genotypes
    tuple val(sample), path("${sample}_targets.csv"), emit: targets
    tuple val(sample), path("genotype.log"), emit: log
 
    script:
    def negative = negative_cells ? "--negative_cells ${negative_cells} --sample ${sample}" : ""
    """
    python ${baseDir}/bin/long_reads/genotyping.py --input ${counts} --output ${sample}.csv --targets_output ${sample}_targets.csv --error_floor ${params.error_floor} --min_mut_umis ${params.min_mut_umis} --snv_min_base_quality ${params.snv_min_base_quality} ${negative}
    """

    stub:
    """
    touch ${sample}.csv
    touch ${sample}_targets.csv
    touch genotype.log
    """
 
}