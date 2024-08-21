// SCMSEQ_NEW, first version
nextflow.enable.dsl = 2

// Include here
include { long_reads } from "./subworkflows/long_reads/main"


//


//----------------------------------------------------------------------------//
// SCM-seq pipeline entry points
//----------------------------------------------------------------------------//

//

workflow LONG_READS {

    // Create a Channel of nanopore fastqs from the input csv file
    ch_fastqs = Channel.fromPath(params.input_fastqs).splitCsv(header: true).map { row -> [row.sample, row.lane, row.fastq]}
    long_reads(ch_fastqs)
    long_reads.out.results.view()

}

//

// Default message
workflow {
    
    println "\n"
    println "Hi there! This is the new version of the SCM-seq pre-processing toolkit. The LONG_READS entry point is currently supported."
    println "Usage: nextflow run main.nf -c <config> -params-file <params> -profile <profile> -entry LONG_READS"
    println "See https://github.com/dimadatascience/scmseq_new/blob/develop/main.nf ./config and ./params for configurations and options available."
    println "N.B. BETA version under active development."
    println "\n"

}
