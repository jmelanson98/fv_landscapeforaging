#!/bin/bash
#SBATCH --job-name=colony_run
#SBATCH --array=0-11
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=logs/colony_%A_%a.out
#SBATCH --error=logs/colony_%A_%a.err

# Load any needed modules
module load gcc

# Your list of .DAT files (replace with real names or glob dynamically)
FILES=(augmented_pp0.DAT augmented_pp0.2.DAT augmented_pp0.4.DAT augmented_pp0.6.DAT augmented_pp0.8.DAT augmented_pp1.DAT
augmented_pp0_poly.DAT augmented_pp0.2_poly.DAT augmented_pp0.4_poly.DAT augmented_pp0.6_poly.DAT augmented_pp0.8_poly.DAT augmented_pp1_poly.DAT)

# Run COLONY on the file for this array task
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}
./colony2s.gnu.out IFN:$INPUT
