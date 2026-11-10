// BLAZE module

nextflow.enable.dsl = 2

//

process BLAZE {

  tag "${sample_name}: ${lane}"
  publishDir "${params.outdir}/${sample_name}/logs/", mode: 'copy', pattern: 'blaze.log'

  input:
  // known_bc: the sample's short-read cell barcodes (plain or gzipped, one per
  // line, Cell Ranger "-1" suffix allowed), or [] to let BLAZE build its own
  // whitelist from the long reads
  tuple val(sample_name), val(lane), path(fastq), path(known_bc)

  output:
  tuple val(sample_name), path("${sample_name}_${lane}_matched_reads.fastq.gz"), emit: filtered_fastq
  tuple val(sample_name), path("blaze.log"), emit: log

  script:
  if (known_bc) {
    // Assign reads to the known barcodes, with no count cutoff. The list is
    // normalised to the bare 16-mers FIXTAGS and the negative-cell join expect.
    """
    zcat -f ${known_bc} | cut -d- -f1 | awk 'NF' | sort -u > known_bc_list.txt

    blaze \
    --no-whitelisting \
    --known-bc-list known_bc_list.txt \
    --threads ${task.cpus} \
    --output-prefix ${sample_name}_${lane}_ \
    --max-edit-distance ${params.blaze_edit_distance} \
    ${fastq}

    echo "Barcode mode: known-bc-list (\$(wc -l < known_bc_list.txt) barcodes)" > blaze.log
    cat .command.out >> blaze.log
    """
  } else {
    """
    blaze \
    --expect-cells ${params.blaze_expect_cells} \
    --threads ${task.cpus} \
    --output-prefix ${sample_name}_${lane}_ \
    --max-edit-distance ${params.blaze_edit_distance} \
    ${fastq}

    cp .command.out blaze.log
    """
  }

  stub:
  def mode = known_bc ? 'known-bc-list' : 'whitelisting'
  """
  touch ${sample_name}_${lane}_matched_reads.fastq.gz
  echo "Barcode mode: ${mode}" > blaze.log
  """

}
