#!/bin/bash
#SBATCH --job-name=nanogen_ladder
#SBATCH --cpus-per-task=2
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --partition=medium
# Attribution ladder S0 -> S6 (+ S0+6) on the PALM14E on-target subset (wayfinder ticket #21).
# Submit from the scripts directory (umi-measurement/, kit-indel/, regression/ side by side):
#   sbatch --output=logs/%x_%j.log regression/run_ladder.sh                 # full run
#   sbatch --partition=short --time=2:00:00 --mem=32gb --output=logs/%x_%j.log \
#          regression/run_ladder.sh --max-cells 600 --out-suffix _smoke       # smoke test
#   sbatch --dependency=afterok:<jobid> --partition=short --mem=4gb --output=logs/%x_%j.log \
#          --wrap "singularity exec -B /hpcnfs <IMG> python regression/fault_review.py $D/regression"  # acceptance 3
# Reads PALM14E_subset in place; writes only to PALM14E_subset/regression<suffix>/ (group-restricted):
# aggregate/ (public-safe, targets as T1..Tn) and private/ (per-cell / per-flip, HPC only).
set -euo pipefail
umask 027
D=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E_subset
P=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E
REF=/hpcnfs/scratch/P_DIMA_SCMSEQ/resources/ref_standard/GRCh38_full_analysis_set_plus_decoy_hla.fa
IMG=/hpcnfs/data/DIMA/singularity_cache/acox1-scmseq-1.0.img
S=$(pwd)

EXTRA=()
SUFFIX=""
while [ $# -gt 0 ]; do
  case "$1" in
    --out-suffix) SUFFIX="$2"; shift 2 ;;
    *) EXTRA+=("$1"); shift ;;
  esac
done
OUT=$D/regression$SUFFIX
mkdir -p "$OUT"

singularity exec -B /hpcnfs $IMG env PYTHONPATH=$S/umi-measurement:$S/kit-indel:$S/regression \
  python $S/regression/ladder.py --cells $D/cells --bed $P/target.bed --ref $REF \
  --cache $D/records.pkl --counts $D/counts_rerun.csv \
  --calls $D/palm14e_rerun.csv $P/results/palm14e/palm14e.csv \
  --out "$OUT" ${EXTRA[@]+"${EXTRA[@]}"}
