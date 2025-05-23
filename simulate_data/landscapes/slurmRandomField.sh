#!/bin/bash
#SBATCH --job-name=landscape_sim
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --array=1-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00

module load StdEnv/2023 gcc r/4.3.1 gdal proj

# List of required packages
packages=("matrixStats" "sp" "gstat" "ggplot2" "reshape2" "raster" "rasterVis" "parallel" "future" "furrr")


# Loop through the list of packages and install them if missing
for package in "${packages[@]}"; do
  Rscript -e "if (!require('$package')) install.packages('$package', repos='https://cloud.r-project.org/')"
done

#Make directory to store outputs
mkdir -p landscapes/random_field_range10

#Run landscape simulation script
Rscript parallelizeSimulateLandscape.R