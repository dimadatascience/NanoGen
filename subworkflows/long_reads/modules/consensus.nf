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
    # Group UMIs per cell x target: primary mapped reads with MAPQ >= min_mapping_quality,
    # one read group per cluster (<target_id>_<cluster ID>, also the MI tag)
    python ${baseDir}/bin/long_reads/umi_grouping.py \
    -i ${bam} \
    -b ${params.bedfile} \
    -o grouped.bam \
    --min_mapq ${params.min_mapping_quality} \
    --max_ld ${params.umi_max_ld} \
    --ratio ${params.umi_ratio} \
    --small_reads ${params.umi_small_reads}

    # One bam-readcount -p pass per target type, concatenated SNV pass first: SNV
    # targets keep the base-quality filter, indel targets have none, as it would
    # drop a deletion on its anchor base. A pass with no targets is skipped.
    awk -F'\\t' '\$5 !~ /^[+-]/' ${params.bedfile} > snv_targets.bed
    awk -F'\\t' '\$5 ~ /^[+-]/' ${params.bedfile} > indel_targets.bed
    : > allcounts.count
    if [ -s snv_targets.bed ]; then
        bam-readcount -f ${params.minimap2_input_genome} --site-list snv_targets.bed \
        -b ${params.snv_min_base_quality} grouped.bam -p >> allcounts.count
    fi
    if [ -s indel_targets.bed ]; then
        bam-readcount -f ${params.minimap2_input_genome} --site-list indel_targets.bed \
        -b 0 grouped.bam -p >> allcounts.count
    fi

    # Consensus labels and sub-n_min clusters per target
    python ${baseDir}/bin/long_reads/count_consensus.py \
    -i allcounts.count \
    -t ${params.bedfile} \
    --cell_barcode ${cell} \
    --min_read ${params.min_reads} \
    --min_fraction ${params.min_fraction}

    rm grouped.bam* snv_targets.bed indel_targets.bed allcounts.count
    """

    stub:
    """
    touch ${cell}_table.csv
    """

}
