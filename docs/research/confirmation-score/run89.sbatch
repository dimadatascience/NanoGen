#!/bin/bash
#SBATCH -J q89_score
#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
# Ticket 89 scoring. Usage: sbatch --chdir=<q89 dir> run89.sbatch <stage1|stage2>
# Python from the pinned v3 image; call() imported from genotype-v3@61d8809. Aggregate output only.
set -euo pipefail
umask 027
Q=/beegfs/scratch/ieo6129/scmseq/q89_score
IMG=/hpcnfs/scratch/DIMA/zhan/scratch/nanogen_img_check/nxf_cache/ghcr.io-dimadatascience-nanogen@sha256-764368d67cdc01c3e40a83157154b97096b7b6ffb2ec740c4c8a4a408720f9ed.img
GT=/beegfs/scratch/ieo6129/scmseq/q65_probe/code/NanoGen-v3/bin/long_reads
echo "== md5"; md5sum $GT/genotyping.py $Q/code/*.py
mkdir -p $Q/out
RUN="singularity exec -B /hpcnfs -B /beegfs --env GENOTYPING_DIR=$GT,PYTHONPATH=$Q/code $IMG python3 -B"
NOBC='[ACGT]{16}'
case $1 in
  stage1) $RUN $Q/code/nosignal89.py $Q/out 2>&1 | grep -v -E "$NOBC" ;;
  stage2) $RUN $Q/code/score89.py $Q/out 2>&1 | grep -v -E "$NOBC" ;;
esac
echo "== done"
