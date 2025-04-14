#!/bin/bash
#SBATCH --job-name=fit_stan
#SBATCH --output=logs/stan_%A_%a.out
#SBATCH --error=logs/stan_%A_%a.err
#SBATCH --array=1-150
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

module load StdEnv/2023 gcc r/4.3.1 gdal proj

# List of required packages
packages=("matrixStats" "sp" "gstat" "ggplot2" "reshape2" "raster" "rasterVis" "parallel" "future" "furrr" "dplyr" "tidyr" "gridExtra")


# Loop through the list of packages and install them if missing
for package in "${packages[@]}"; do
  Rscript -e "if (!require('$package')) install.packages('$package', repos='https://cloud.r-project.org/')"
done

Rscript testPopeModel.R