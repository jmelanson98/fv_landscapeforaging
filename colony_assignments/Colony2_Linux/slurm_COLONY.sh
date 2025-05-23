#!/bin/bash
#SBATCH --job-name=colony_run
#SBATCH --array=0-4         # Adjust for number of .DAT files
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G            # Adjust depending on dataset size
#SBATCH --time=24:00:00
#SBATCH --output=colony_%A_%a.out
#SBATCH --error=colony_%A_%a.err

# Load any needed modules
module load gcc

# Your list of .DAT files (replace with real names or glob dynamically)
FILES=(impatiens2022.DAT impatiens2022_1.DAT impatiens2023_1.DAT mixtus2022_1.DAT mixtus2023_1.DAT)

# Run COLONY on the file for this array task
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}
./colony2s.gnu.out IFN:$INPUT
