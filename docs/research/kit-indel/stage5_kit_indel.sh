#!/bin/bash
#SBATCH --job-name=nanogen_s5
#SBATCH --cpus-per-task=2
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=stage5.log
#SBATCH --partition=medium
# KIT indel labelling: main's rule vs the decided tolerant rule (k = 0, 1, 2), near-miss spectrum and
# the site-line blind spot, on main's grouping (stage1 allcounts.count) and the new grouping
# (records.pkl, hybrid T=3) -> kit_out/.
D=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E_subset
P=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E
REF=/hpcnfs/scratch/P_DIMA_SCMSEQ/resources/ref_standard/GRCh38_full_analysis_set_plus_decoy_hla.fa
cd $D
singularity exec -B /hpcnfs /hpcnfs/data/DIMA/singularity_cache/acox1-scmseq-1.0.img env PYTHONPATH=scripts \
  python scripts/kit_indel.py --cells cells --bed $P/target.bed --ref $REF --counts counts_rerun.csv --calls palm14e_rerun.csv \
  --cache records.pkl --out kit_out
