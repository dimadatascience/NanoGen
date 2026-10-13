#!/bin/bash
#SBATCH --job-name=gate_mock_bam
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=2:00:00
# Merge gate (#41) mocks: the tagged on-target PALM14E BAM restricted to N random cells (seed 41),
# to size the v2 run and the gate. Usage: sbatch make_mock.sh <gate dir> <N> [<N> ...]
# Writes <gate dir>/mock<N>/palm14e_mock<N>.bam and samplesheet.csv (group-restricted).
set -euo pipefail
umask 007
G=$1; shift
D=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E_subset
IMG=/hpcnfs/scratch/DIMA/zhan/scratch/nanogen_img_check/nxf_cache/ghcr.io-dimadatascience-nanogen@sha256-764368d67cdc01c3e40a83157154b97096b7b6ffb2ec740c4c8a4a408720f9ed.img
awk -F, 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "cell") c = i; if (!c) exit 1; next} {print $c}' \
  "$D/palm14e_rerun.csv" | sort -u > "$G/all_cells.txt"
echo "cells available: $(wc -l < "$G/all_cells.txt")"
for N in "$@"; do
  M=$G/mock$N
  mkdir -p "$M"
  shuf -n "$N" --random-source=<(yes 41) "$G/all_cells.txt" > "$M/cells.txt"
  singularity exec -B /hpcnfs "$IMG" samtools view -@ 2 -b -D CB:"$M/cells.txt" \
    -o "$M/palm14e_mock$N.bam" "$D/ontarget/palm14e_ontarget.bam"
  printf 'sample,bam\npalm14e,%s\n' "$M/palm14e_mock$N.bam" > "$M/samplesheet.csv"
  echo "mock$N records: $(singularity exec -B /hpcnfs "$IMG" samtools view -c "$M/palm14e_mock$N.bam")"
done
