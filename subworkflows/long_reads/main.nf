// long_reads subworkflow

// Include here
nextflow.enable.dsl = 2
include { BLAZE } from "./modules/blaze.nf"
include { MERGE_READS } from "./modules/merge_reads.nf"
include { MINIMAP2 } from "./modules/minimap2.nf"
include { FIXTAGS } from "./modules/fix_bam_tags.nf"
include { SPLIT_BAM } from "./modules/split_bam.nf"
include { CONSENSUS } from "./modules/consensus.nf"
include { COLLAPSE } from "./modules/collapse_consensus.nf"
include { GENOTYPING } from "./modules/genotyping.nf"


// 

process publish_long_reads {

    publishDir "${params.outdir}/${sample_name}/", mode: 'copy'

    input:
    tuple val(sample_name), ...

    output:
    ...

    script:
    """
    echo moving everything to ${params.sc_outdir}
    """

}
 
// 

//----------------------------------------------------------------------------//
// long_reads subworkflow
//----------------------------------------------------------------------------//

workflow long_reads {
     
    take:
        ch_input 

    main:

        BLAZE(ch_input)
        MERGE_READS(BLAZE.out.filtered_fastq.groupTuple(by:0))
        MINIMAP2(MERGE_READS.out.filtered_fastq)
        FIXTAGS(MINIMAP2.out.bam)
        SPLIT_BAM(FIXTAGS.out.bam)
        ch_cell_bams = SPLIT_BAM.out.cell_bams
            .map { it ->
                def sample = it[0] 
                def paths = it[1]     
                return paths.collect { cell_path ->
                    def path_splitted = cell_path.toString().split('/')
                    def cell = path_splitted[-1].toString().split('\\.')[0]
                    return [sample, cell, cell_path]
                }
            } 
            .flatMap { it } 
        CONSENSUS(ch_cell_bams)
        COLLAPSE(CONSENSUS.out.tables.groupTuple(by: 0))
        GENOTYPING(COLLAPSE.out.allele_counts)

    emit:

        results = GENOTYPING.out.genotypes

}

//----------------------------------------------------------------------------//