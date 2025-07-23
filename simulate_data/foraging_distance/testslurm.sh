#!/bin/bash
##SBATCH --job-name=test
#SBATCH --output=logs/widelimits_%A_%a.out
#SBATCH --error=logs/widelimits_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

module load StdEnv/2023 gcc r/4.3.1 gdal proj

# Define your personal R library path
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.3


Rscript 8_observed_vs_unobserved.R
