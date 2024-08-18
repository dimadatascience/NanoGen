// long_reads subworkflow

// Include here
nextflow.enable.dsl = 2
include { BLAZE } from "./modules/blaze.nf"
include { MERGE_READS } from "./modules/merge_reads.nf"
include { MINIMAP2 } from "./modules/minimap2.nf"
include { FIXTAGS } from "./modules/fix_bam_tags.nf"
include { SPLIT } from "./modules/filter_bam.nf"
// include { MERGE_BAM } from "./modules/merge_bams.nf"
// include { EXTRACT_FASTA } from "../common/modules/extract_fasta.nf"
// include { SPLIT_BARCODES } from "../common/modules/split_barcodes.nf"


// 

// process publish_long_reads {
// 
//     publishDir "${params.outdir}/${sample_name}/", mode: 'copy'
// 
//     input:
//     tuple val(sample_name), ...
// 
//     output:
//     ...
// 
//     script:
//     """
//     echo moving everything to ${params.sc_outdir}
//     """
// 
// }
 
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

    emit:

        results = FIXTAGS.out.bam
        // results = MERGE_R1.out.R1

}

//----------------------------------------------------------------------------//