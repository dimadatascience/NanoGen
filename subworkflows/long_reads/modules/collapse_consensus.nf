// COLLAPSE  module
 
nextflow.enable.dsl = 2
 
//


process COLLAPSE {
 
    tag "${sample}"
 
    input:
        tuple val(sample), path(files)
 
    output:
        tuple val(sample), path("counts.tsv.gz"), emit: allele_counts
 
    script:
    """
    outfile="counts.tsv"
    files=(${files})
    cat "\${files[0]}" > \$outfile
    for f in "\${files[@]:1}"; do
        tail -n +2 "\$f" >> \$outfile
    done
    gzip --fast \$outfile
    """

    stub:
    """
    touch counts.tsv.gz
    """
 
}