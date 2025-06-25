#!/bin/bash

#SBATCH --job-name=hvgs

#SBATCH --cpus-per-task 1

#SBATCH --mem 300G

#SBATCH --partition=long

#SBATCH --output=submit_hvgs.txt

#SBATCH --error=submit_hvgs.err

module load R/4.4.2-foss-2021b

Rscript 03_HVGs_snRNA.R
