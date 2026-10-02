#!/bin/bash
#SBATCH --job-name=nanogen_s1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=12:00:00
#SBATCH --output=stage1.log
#SBATCH --partition=medium
# Re-run main@c442a32 per-cell consensus on the on-target subset with the
# baseline params (nextflow.config defaults), keeping intermediates.
set -euo pipefail
D=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E_subset
P=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E
REF=/hpcnfs/scratch/P_DIMA_SCMSEQ/resources/ref_standard/GRCh38_full_analysis_set_plus_decoy_hla.fa
IMG=/hpcnfs/data/DIMA/singularity_cache/acox1-scmseq-1.0.img
export D P REF IMG
cd $D
mkdir -p cells split
run() { singularity exec -B /hpcnfs $IMG "$@"; }
export -f run

# split by CB exactly like SPLIT_BAM
if [ ! -f split/.done ]; then
  cd split
  run samtools split -M -1 -@ 16 -d CB -f '%!.bam' $D/ontarget/palm14e_ontarget.bam
  touch .done
  cd ..
fi

one() {
  bam=$1; cell=$(basename $bam .bam)
  mkdir -p $D/cells/$cell && cd $D/cells/$cell
  run samtools sort -o $cell.bam $bam 2>/dev/null
  run bash $D/NanoGen/bin/long_reads/make_consensus.sh -i $cell.bam -b $P/target.bed -r $REF \
      -q 20 -c 2 -w 100 -m 20 -d 1 > mc.log 2>&1 || { echo "FAIL $cell" ; return 0; }
  run python $D/NanoGen/bin/long_reads/count_consensus.py -i temp_$cell/allcounts.count -t $P/target.bed \
      --cell_barcode $cell --min_read 2 --min_fraction 0.6 >> mc.log 2>&1 || echo "FAILCC $cell"
  # keep: filtered_input.bam (BX/UG tagged), grouped_reads.tsv, allcounts.count, table
  mv temp_$cell/allcounts.count . 2>/dev/null || true
  rm -rf temp_$cell output.bam output.sam bam_header.sam $cell.bam
}
export -f one
ls $D/split/*.bam | xargs -P 16 -I{} bash -c 'one {}'

# collapse + genotype with main's genotyping.py
first=1
for f in cells/*/*_table.csv; do
  if [ $first = 1 ]; then cat $f; first=0; else tail -n +2 $f; fi
done > counts_rerun.csv
run python $D/NanoGen/bin/long_reads/genotyping.py --input counts_rerun.csv --output palm14e_rerun.csv --baseline_error 0.0001
echo "cells processed: $(ls -d cells/* | wc -l)"
echo DONE
