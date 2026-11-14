#!/bin/bash
#SBATCH -J q94_v3ub
#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
# Ticket 94: per-cell v3-ub calls. Usage: sbatch --chdir=<q94 dir> run94.sbatch
# Library = q89_score/code @ 23e407c; call() from genotype-v3@61d8809; inputs q79_full/runs/v3_<S>. Aggregate log only.
set -euo pipefail
umask 027
Q=/beegfs/scratch/ieo6129/scmseq/q94_v3ub
L=/beegfs/scratch/ieo6129/scmseq/q89_score
IMG=/hpcnfs/scratch/DIMA/zhan/scratch/nanogen_img_check/nxf_cache/ghcr.io-dimadatascience-nanogen@sha256-764368d67cdc01c3e40a83157154b97096b7b6ffb2ec740c4c8a4a408720f9ed.img
GT=/beegfs/scratch/ieo6129/scmseq/q65_probe/code/NanoGen-v3/bin/long_reads
echo "== md5"; md5sum $GT/genotyping.py $L/code/*.py $Q/code/*.py
mkdir -p $Q/out
singularity exec -B /hpcnfs -B /beegfs --env GENOTYPING_DIR=$GT,PYTHONPATH=$L/code $IMG \
  python3 -B $Q/code/mat94.py $Q/out $L/out 2>&1 | grep -v -E '[ACGT]{16}'
echo "== done"
