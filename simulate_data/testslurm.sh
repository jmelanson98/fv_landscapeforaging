#!/bin/bash
##SBATCH --job-name=test
#SBATCH --output=logs/traplevel_%A_%a.out
#SBATCH --error=logs/traplevel_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=12G
#SBATCH --time=08:00:00

module load StdEnv/2023 gcc r/4.3.1 gdal proj

# Define your personal R library path
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.3


Rscript 10_trap_landscape_effects.R
