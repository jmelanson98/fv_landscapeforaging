#!/bin/bash
#SBATCH --job-name=colony_run
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=logs/colony_%A_%a.out
#SBATCH --error=logs/colony_%A_%a.err

# Load any needed modules
module load gcc

# Your list of .DAT files (replace with real names or glob dynamically)
FILES=(impatiens_2023wq.DAT mixtus_2023wq.DAT impatiens_2022.DAT impatiens_2023.DAT mixtus_2022.DAT mixtus_2023.DAT)

# Run COLONY on the file for this array task
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}
./colony2s.gnu.out IFN:$INPUT
