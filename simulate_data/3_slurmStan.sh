#!/bin/bash
#SBATCH --job-name=fit_stan
#SBATCH --output=logs/stanGQ_%A_%a.out
#SBATCH --error=logs/stanGQ_%A_%a.err
#SBATCH --array=1-150
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00

module load StdEnv/2023 gcc r/4.3.1 gdal proj

# List of required packages
packages=("matrixStats" "sp" "gstat" "ggplot2" "reshape2" "raster" "rasterVis" "parallel" "future" "furrr" "dplyr" "tidyr" "gridExtra")

# Define your personal R library path
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.3

# Create the personal R library directory if it doesn't exist
mkdir -p "$R_LIBS_USER"

# Loop through the list of packages and install them if missing
for package in "${packages[@]}"; do
  Rscript -e "if (!require('$package')) install.packages('$package', repos='https://cloud.r-project.org/')"
done

Rscript 3A_testGeneratedQuantitiesModel.R
