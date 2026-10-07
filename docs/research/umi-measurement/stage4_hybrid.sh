#!/bin/bash
#SBATCH --job-name=nanogen_s4
#SBATCH --cpus-per-task=2
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=stage4.log
#SBATCH --partition=medium
# Hybrid UMI merge rule (n_b <= T clause) vs plain directional LD<=2, from records.pkl -> hybrid_out/.
D=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E_subset
P=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E
cd $D
singularity exec -B /hpcnfs /hpcnfs/data/DIMA/singularity_cache/acox1-scmseq-1.0.img env PYTHONPATH=scripts \
  python scripts/hybrid.py --bed $P/target.bed --cache records.pkl --out hybrid_out
