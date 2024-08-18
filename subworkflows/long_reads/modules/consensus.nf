// CONSENSUS module

nextflow.enable.dsl = 2

//

process CONSENSUS {

    tag "${sample_name}: ${cell}"

    input:
    tuple val(sample_name), val(cell), path(bam)

    output:
    tuple val(sample_name), path("${cell}_table.csv"), emit: tables

    script:
    """ 
    bash ${baseDir}/bin/long_reads/split_bam_by_tag_and_window_count_consensus.sh \
    -i ${bam} \
    -b ${params.bedfile} \
    -r ${params.minimap2_input_genome} \
    -q ${params.min_mapping_quality} \
    -c ${params.min_reads} \
    -w ${params.window_width} \
    -m ${params.min_base_quality} \
    -d ${params.min_dist_umi_collapsing}

    # Problems with realtive script paths in split_bam_by_tag_and_window_count_consensus: FIX IT
    """

    stub:
    """
    touch ${cell}_table.csv
    """

}
