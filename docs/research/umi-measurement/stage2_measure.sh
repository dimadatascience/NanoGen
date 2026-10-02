#!/bin/bash
#SBATCH --job-name=nanogen_s2
#SBATCH --cpus-per-task=2
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=stage2.log
#SBATCH --partition=medium
# Run the UMI measurements on stage-1 intermediates; writes aggregate CSVs to measure_out/.
D=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E_subset
P=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E
cd $D
singularity exec -B /hpcnfs /hpcnfs/data/DIMA/singularity_cache/acox1-scmseq-1.0.img \
  python scripts/measure_umi.py --cells cells --bed $P/target.bed --out measure_out
