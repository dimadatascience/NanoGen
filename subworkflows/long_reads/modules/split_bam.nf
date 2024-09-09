// SPLIT_BAM module

nextflow.enable.dsl = 2

//

process SPLIT_BAM {

    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}/logs/", mode: 'copy', pattern: 'target.log'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path('*.bam'), emit: cell_bams
    tuple val(sample_name), path('target.log'), emit: log

    script:
    """ 
    # Index
    samtools sort -@ ${task.cpus} ${bam} -o sorted_${bam}
    samtools index sorted_${bam}

    # Count total reads (both mapped and unmapped)
    total_reads=\$(samtools view sorted_${bam} | wc -l)
    # Count reads on target regions using bedtools intersect
    reads_on_target=\$(bedtools intersect -a sorted_${bam} -b ${params.bedfile} -u |  samtools view - | wc -l)

    # Get statistics
    echo "On target: \$reads_on_target" > target.log
    echo "Total: \$total_reads" >> target.log
    bedtools multicov -bams sorted_${bam} -bed ${params.bedfile} >> target.log 

    samtools split -M -1 -@ ${task.cpus} -d CB -f '%!.bam' ${bam}
    rm -f sorted_${bam} ${bam}
    
    for bam_file in *.bam; do
        samtools sort -@ ${task.cpus} -o sorted_\${bam_file} \${bam_file}
        mv sorted_\${bam_file} \${bam_file}
    done

    """

    stub:
    """
    python ${baseDir}/bin/long_reads/random_mock.py \
    touch target.log
    rm -f ${bam}
    """

}
