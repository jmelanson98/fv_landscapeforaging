#!/bin/bash
##SBATCH --job-name=test
#SBATCH --output=methods_comparison/landscape_effects/logs/landscapesensitivity_%A_%a.out
#SBATCH --error=methods_comparison/landscape_effects/logs/landscapesensitivity_%A_%a.err
#SBATCH --array=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00

module load StdEnv/2023 gcc r/4.3.1 gdal proj

# Define your personal R library path
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.3


Rscript 13_landscape_sensitivity_fit.R
