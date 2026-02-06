#!/bin/bash
##SBATCH --job-name=makeraster
#SBATCH --output=logs/makeraster_%A_%a.out
#SBATCH --error=logs/makeraster_%A_%a.err
#SBATCH --array=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=12:00:00

module load StdEnv/2023 gcc r/4.3.1 gdal proj

# Define your personal R library path
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.3


Rscript createraster.R
