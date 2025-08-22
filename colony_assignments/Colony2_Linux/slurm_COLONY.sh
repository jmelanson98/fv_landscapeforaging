#!/bin/bash
#SBATCH --job-name=colony_run
#SBATCH --array=10-14
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=logs/colony_%A_%a.out
#SBATCH --error=logs/colony_%A_%a.err

# Load any needed modules
module load gcc

# Your list of .DAT files (replace with real names or glob dynamically)
FILES=(impatiens2022_final1.DAT impatiens2022_final2.DAT impatiens2022_final3.DAT impatiens2022_final4.DAT impatiens2022_final5.DAT
impatiens2023_final1.DAT impatiens2023_final2.DAT impatiens2023_final3.DAT impatiens2023_final4.DAT impatiens2023_final5.DAT
mixtus2022_final1.DAT mixtus2022_final2.DAT mixtus2022_final3.DAT mixtus2022_final4.DAT mixtus2022_final5.DAT
mixtus2023_final1.DAT mixtus2023_final2.DAT mixtus2023_final3.DAT mixtus2023_final4.DAT mixtus2023_final5.DAT)

# Run COLONY on the file for this array task
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}
./colony2s.gnu.out IFN:$INPUT
