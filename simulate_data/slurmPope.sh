#!/bin/bash
#SBATCH --job-name=bee_sim
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --array=1-150
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --account=def-youraccount  # Replace with your Alliance account
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END,FAIL

module load r/4.3.1  # Adjust version as needed

Rscript simulate_one_run.R

files <- list.files("results", pattern = "sim_result_.*rds", full.names = TRUE)
results <- lapply(files, readRDS)
param_grid <- readRDS("param_grid.rds")

full_results <- bind_cols(param_grid, tibble(result = results))
