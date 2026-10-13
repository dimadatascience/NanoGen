#!/bin/bash
#SBATCH --job-name=gate_v2_head
#SBATCH --partition=medium
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=4-00:00:00
# Merge gate (#41): genotype-v2 from the tagged on-target PALM14E BAM (--input_type tagged_bam).
# The Nextflow head runs here; processes go to slurm (queue and resources from the environment,
# see merge_gate.config). Usage: sbatch --chdir=<run dir> run_v2.sh <NanoGen checkout> <samplesheet>
# Writes <run dir>/{results,work,trace.txt}: patient data, group-restricted, never committed.
# The Nextflow log holds cell barcodes in process tags: it stays in the run dir, not in this log.
set -euo pipefail
umask 007
NG=$1; SS=$2
P=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E
export JAVA_HOME=/hpcnfs/data/PGP/zhan/miniconda3/envs/biomodal
export JAVA_CMD=$JAVA_HOME/bin/java
export PATH=$JAVA_HOME/bin:$PATH
export NXF_SINGULARITY_CACHEDIR=/hpcnfs/scratch/DIMA/zhan/scratch/nanogen_img_check/nxf_cache
export NXF_OPTS="-Xms1g -Xmx6g"
echo "host: $(hostname)  queue: ${GATE_QUEUE:-short}  v2: $(cat "$NG/.gate_rev")"
status=0
nextflow run "$NG/main.nf" -entry LONG_READS -profile singularity \
  -c "$P/nextflow.config" -c "$NG/docs/research/regression/merge_gate.config" \
  --input "$SS" --input_type tagged_bam --outdir "$PWD/results" \
  -work-dir "$PWD/work" -ansi-log false -resume > nextflow_head.log 2>&1 || status=$?
echo "--- process status"
awk -F'\t' 'NR>1{split($3,a," "); print a[1], $4}' trace.txt | sort | uniq -c
echo "nextflow exit: $status"
exit $status
