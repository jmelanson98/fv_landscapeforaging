#!/bin/bash
##SBATCH --job-name=floralforaging
#SBATCH --output=logs/floralforaging_%A_%a.out
#SBATCH --error=logs/floralforaging_%A_%a.err
#SBATCH --array=1-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=72:00:00

module load StdEnv/2023 gcc r/4.3.1 gdal proj

# Define your personal R library path
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.3


Rscript floral_foraging.R
