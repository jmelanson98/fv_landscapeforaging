#!/bin/bash
#SBATCH --job-name=fit_cmdstanr
#SBATCH --output=logs/cmdstanr_%A_%a.out
#SBATCH --error=logs/cmdstanr_%A_%a.err
#SBATCH --array=74
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=12G
#SBATCH --time=10:00:00

module load StdEnv/2023 gcc r/4.3.1 gdal proj

# Define your personal R library path
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.3

Rscript 8A_reducesum_unobserved.R

