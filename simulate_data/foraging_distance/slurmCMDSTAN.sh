#!/bin/bash
#SBATCH --job-name=fit_cmdstanr
#SBATCH --output=logs/trapobserved_landscape_%A_%a.out
#SBATCH --error=logs/trapobserved_landscape_%A_%a.err
#SBATCH --array=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=12G
#SBATCH --time=04:00:00

module load StdEnv/2023 gcc r/4.3.1 gdal proj

# Define your personal R library path
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.3
export STAN_NUM_THREADS=4

Rscript 10_trap_landscape_effects.R

