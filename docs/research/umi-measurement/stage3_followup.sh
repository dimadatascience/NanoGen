#!/bin/bash
#SBATCH --job-name=nanogen_s3
#SBATCH --cpus-per-task=2
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=stage3.log
#SBATCH --partition=medium
# Follow-up checks on merge validity; caches per-read table as records.pkl (HPC only).
D=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E_subset
P=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E
cd $D
singularity exec -B /hpcnfs /hpcnfs/data/DIMA/singularity_cache/acox1-scmseq-1.0.img env PYTHONPATH=scripts \
  python scripts/followup.py --cells cells --bed $P/target.bed --out measure_out --cache records.pkl
