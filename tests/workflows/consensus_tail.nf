// Test-only workflow: the CONSENSUS -> COLLAPSE -> GENOTYPING tail of
// long_reads, fed per-cell SAM text converted to BAM. Once the fixture
// generator (#27) lands, point the test at its BAMs and drop SAM_TO_BAM.

nextflow.enable.dsl = 2

include { CONSENSUS }  from "../../subworkflows/long_reads/modules/consensus.nf"
include { COLLAPSE }   from "../../subworkflows/long_reads/modules/collapse_consensus.nf"
include { GENOTYPING } from "../../subworkflows/long_reads/modules/genotyping.nf"

process SAM_TO_BAM {

    tag "${sample}: ${cell}"

    input:
    tuple val(sample), val(cell), path(sam)

    output:
    tuple val(sample), val(cell), path("${cell}.bam")

    script:
    """
    samtools view -b -o ${cell}.bam ${sam}
    """

}

workflow CONSENSUS_TAIL {

    take:
        ch_sams // [sample, cell, sam]

    main:
        SAM_TO_BAM(ch_sams)
        CONSENSUS(SAM_TO_BAM.out)
        // Sort cell tables by name so the collapsed table is reproducible
        COLLAPSE(CONSENSUS.out.tables.groupTuple(by: 0, sort: { a, b -> a.name <=> b.name }))
        GENOTYPING(COLLAPSE.out.allele_counts)

    emit:
        counts    = COLLAPSE.out.allele_counts
        genotypes = GENOTYPING.out.genotypes
        gt_log    = GENOTYPING.out.log

}
