#!/bin/bash
#SBATCH -J pool_guard
#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 00:20:00
#SBATCH -o pool_guard.%j.log
cd /beegfs/scratch/ieo6129/scmseq/q83_pool_guard/
singularity exec -B /hpcnfs /hpcnfs/data/DIMA/singularity_cache/acox1-scmseq-1.0.img python3 pool_guard.py eb.tsv 100000
