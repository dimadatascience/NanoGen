// SCMSEQ_NEW, first commit
nextflow.enable.dsl = 2

// Include here
include { long_reads } from "./subworkflows/long_reads/main"


//


// Create a Channel of nanopore fastqs from the input csv file
ch_fastqs = Channel.fromPath(params.input_fastqs).splitCsv(header: true).map { row -> [row.sample, row.lane, row.fastq]}


//


//----------------------------------------------------------------------------//
// SCM-seq pipeline entry points
//----------------------------------------------------------------------------//

//

workflow LONG_READS {

    long_reads(ch_fastqs)
    long_reads.out.results.view()

}

//

// Mock
workflow {
    
    Channel.of(1,2,3,4) | view

}
