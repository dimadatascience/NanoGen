#!/bin/bash
#SBATCH -J q91_indel
#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
# Ticket 91 indel arm. Usage: sbatch --chdir=<q91 dir> run91.sbatch <score|diag>
# Python from the pinned v3 image; call() imported from genotype-v3@61d8809; library = q89 code @ 23e407c.
set -euo pipefail
umask 027
Q=/beegfs/scratch/ieo6129/scmseq/q91_indel
Q89=/beegfs/scratch/ieo6129/scmseq/q89_score
IMG=/hpcnfs/scratch/DIMA/zhan/scratch/nanogen_img_check/nxf_cache/ghcr.io-dimadatascience-nanogen@sha256-764368d67cdc01c3e40a83157154b97096b7b6ffb2ec740c4c8a4a408720f9ed.img
GT=/beegfs/scratch/ieo6129/scmseq/q65_probe/code/NanoGen-v3/bin/long_reads
echo "== md5"; md5sum $GT/genotyping.py $Q89/code/*.py $Q/code/*.py
mkdir -p $Q/out
RUN="singularity exec -B /hpcnfs -B /beegfs --env GENOTYPING_DIR=$GT,PYTHONPATH=$Q/code:$Q89/code $IMG python3 -B"
NOBC='[ACGT]{16}'
case $1 in
  score) $RUN $Q/code/score91.py $Q89/out $Q/out 2>&1 | grep -v -E "$NOBC" ;;
  diag)  $RUN $Q/code/diag91.py $Q/out 2>&1 | grep -v -E "$NOBC" ;;
esac
echo "== done"
