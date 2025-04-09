#!/bin/bash
#SBATCH --job-name=bee_sim
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --array=1-150
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00

module load r/4.3.1  # Adjust version as needed

Rscript parallelizePopeSim.R


