// Test-only workflow: the CONSENSUS -> COLLAPSE -> GENOTYPING tail of
// long_reads, fed per-cell BAMs as SPLIT_BAM emits them.

nextflow.enable.dsl = 2

include { CONSENSUS }  from "../../subworkflows/long_reads/modules/consensus.nf"
include { COLLAPSE }   from "../../subworkflows/long_reads/modules/collapse_consensus.nf"
include { GENOTYPING } from "../../subworkflows/long_reads/modules/genotyping.nf"

workflow CONSENSUS_TAIL {

    take:
        ch_cell_bams // [sample, cell, bam]

    main:
        CONSENSUS(ch_cell_bams)
        // Sort cell tables by name so the collapsed table is reproducible
        COLLAPSE(CONSENSUS.out.tables.groupTuple(by: 0, sort: { a, b -> a.name <=> b.name }))
        GENOTYPING(COLLAPSE.out.allele_counts)

    emit:
        counts    = COLLAPSE.out.allele_counts
        genotypes = GENOTYPING.out.genotypes
        targets   = GENOTYPING.out.targets
        gt_log    = GENOTYPING.out.log

}
