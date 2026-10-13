#!/bin/bash
#SBATCH --job-name=gate_compare
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
# Merge gate (#41): G0 tie, a1 grouping, a2 caller, a3 labels, G4 migration numbers.
# Usage: sbatch --chdir=<run dir> run_merge_gate.sh <NanoGen checkout> [--max-cells N]
# Reads the v2 run in <run dir> (results/, SPLIT_BAM's per-cell BAMs in work/); writes <run dir>/gate.
set -euo pipefail
umask 007
NG=$1; shift
D=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E_subset
P=$D/../PALM14E
REF=/hpcnfs/scratch/P_DIMA_SCMSEQ/resources/ref_standard/GRCh38_full_analysis_set_plus_decoy_hla.fa
IMG=/hpcnfs/scratch/DIMA/zhan/scratch/nanogen_img_check/nxf_cache/ghcr.io-dimadatascience-nanogen@sha256-764368d67cdc01c3e40a83157154b97096b7b6ffb2ec740c4c8a4a408720f9ed.img
H=$(awk -F'\t' '$3 ~ /:SPLIT_BAM / && $4 == "COMPLETED" {print $2}' trace.txt | tail -1)
[ -n "$H" ] || { echo "no completed SPLIT_BAM in trace.txt"; exit 1; }
CELLS=$(ls -d work/$H*)
[ "$(ls "$CELLS"/*.bam | wc -l)" -gt 0 ] || { echo "no per-cell BAMs"; exit 1; }
echo "v2: $(cat "$NG/.gate_rev")  per-cell BAMs: $(ls "$CELLS"/*.bam | wc -l)"
R=$NG/docs/research
export TMPDIR=${SLURM_TMPDIR:-/tmp}
singularity exec -B /hpcnfs "$IMG" env PYTHONPATH=$R/umi-measurement:$R/kit-indel:$R/regression \
  python "$R/regression/merge_gate.py" --cell-bams "$CELLS" --v2-bin "$NG/bin/long_reads" \
  --v2-results results/palm14e --bed "$P/target.bed" --ref "$REF" \
  --main-calls "$D/palm14e_rerun.csv" --out gate --threads "$SLURM_CPUS_PER_TASK" "$@"
