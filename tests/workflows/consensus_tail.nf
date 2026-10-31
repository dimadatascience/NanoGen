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
        // Negative-control cells as the subworkflow passes them: the file, or [] for none
        ch_negative_cells = params.negative_cells
            ? file(params.negative_cells, checkIfExists: true)
            : []
        GENOTYPING(COLLAPSE.out.allele_counts, ch_negative_cells)

    emit:
        counts    = COLLAPSE.out.allele_counts
        genotypes = GENOTYPING.out.genotypes
        targets   = GENOTYPING.out.targets
        gt_log    = GENOTYPING.out.log

}
