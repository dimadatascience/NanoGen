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
    # Generate cell-target-site-UMI consensus allele counts. The window (-w) and
    # UMI distance (-d) are fixed at their former defaults until cell x target
    # grouping replaces umi_tools and the windows.
    bash ${baseDir}/bin/long_reads/make_consensus.sh \
    -i ${bam} \
    -b ${params.bedfile} \
    -r ${params.minimap2_input_genome} \
    -q ${params.min_mapping_quality} \
    -c ${params.min_reads} \
    -w 100 \
    -m ${params.snv_min_base_quality} \
    -d 1


    # Extract counts from bam-readcount output!!!!
    python ${baseDir}/bin/long_reads/count_consensus.py \
    -i temp_${cell}/allcounts.count \
    -t ${params.bedfile} \
    --cell_barcode ${cell} \
    --min_read ${params.min_reads} \
    --min_fraction ${params.min_fraction}
    

    # Clean up
    rm bam_header.sam
    rm filtered_input.bam*
    rm grouped_reads.tsv
    rm output.bam
    rm output.sam
    rm -r temp_${cell}
    """

    stub:
    """
    touch ${cell}_table.csv
    """

}
