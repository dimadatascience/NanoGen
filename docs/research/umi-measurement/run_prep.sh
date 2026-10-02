#!/bin/bash
#SBATCH --job-name=nanogen_prep
#SBATCH --ntasks=2
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH --output=prep.log
#SBATCH --partition=medium
source ~/.bashrc
conda activate dima
nextflow run NanoGen/prep.nf -c ../PALM14E/nextflow.config -profile singularity \
  --input ../PALM14E/rerun.csv --outdir ontarget -w work -resume
