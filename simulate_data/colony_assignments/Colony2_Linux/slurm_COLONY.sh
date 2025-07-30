#!/bin/bash
#SBATCH --job-name=colony_run
#SBATCH --array=0-47
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/colony_%A_%a.out
#SBATCH --error=logs/colony_%A_%a.err

# Load any needed modules
module load gcc

# Your list of .DAT files (replace with real names or glob dynamically)
FILES=(mixtus_set3_sub1.DAT mixtus_set3_sub0.8.DAT mixtus_set3_sub0.6.DAT mixtus_set3_sub0.4.DAT
mixtus_set3_sub0.2.DAT mixtus_set3_sub0.1.DAT mixtus_set3_sub0.05.DAT mixtus_set3_sub0.025.DAT
mixtus_set2_sub1.DAT mixtus_set2_sub0.8.DAT mixtus_set2_sub0.6.DAT mixtus_set2_sub0.4.DAT
mixtus_set2_sub0.2.DAT mixtus_set2_sub0.1.DAT mixtus_set2_sub0.05.DAT mixtus_set2_sub0.025.DAT
mixtus_set1_sub1.DAT mixtus_set1_sub0.8.DAT mixtus_set1_sub0.6.DAT mixtus_set1_sub0.4.DAT
mixtus_set1_sub0.2.DAT mixtus_set1_sub0.1.DAT mixtus_set1_sub0.05.DAT mixtus_set1_sub0.025.DAT
impatiens_set3_sub1.DAT impatiens_set3_sub0.8.DAT impatiens_set3_sub0.6.DAT impatiens_set3_sub0.4.DATimpatiens_set3_sub0.2.DAT impatiens_set3_sub0.1.DAT impatiens_set3_sub0.05.DAT impatiens_set3_sub0.025.DATimpatiens_set2_sub1.DAT impatiens_set2_sub0.8.DAT impatiens_set2_sub0.6.DAT impatiens_set2_sub0.4.DATimpatiens_set2_sub0.2.DAT impatiens_set2_sub0.1.DAT impatiens_set2_sub0.05.DAT impatiens_set2_sub0.025.DATimpatiens_set1_sub1.DAT impatiens_set1_sub0.8.DAT impatiens_set1_sub0.6.DAT impatiens_set1_sub0.4.DATimpatiens_set1_sub0.2.DAT impatiens_set1_sub0.1.DAT impatiens_set1_sub0.05.DAT impatiens_set1_sub0.025.DAT
)

# Run COLONY on the file for this array task
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}
./colony2s.gnu.out IFN:$INPUT
