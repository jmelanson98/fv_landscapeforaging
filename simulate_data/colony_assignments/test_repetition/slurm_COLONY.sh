#!/bin/bash
#SBATCH --job-name=colony_run
#SBATCH --array=0-49
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=18:00:00
#SBATCH --output=logs/colony_%A_%a.out
#SBATCH --error=logs/colony_%A_%a.err

# Load any needed modules
module load gcc

# Your list of _no_exclusion.DAT files (replace with real names or glob dynamically)
FILES=(mixtus_set1_sub0.6_run1.DAT mixtus_set1_sub0.6_run2.DAT mixtus_set1_sub0.6_run3.DAT mixtus_set1_sub0.6_run4.DAT mixtus_set1_sub0.6_run5.DAT mixtus_set2_sub0.6_run1.DAT mixtus_set2_sub0.6_run2.DAT mixtus_set2_sub0.6_run3.DAT mixtus_set2_sub0.6_run4.DAT mixtus_set2_sub0.6_run5.DAT mixtus_set3_sub0.6_run1.DAT mixtus_set3_sub0.6_run2.DAT mixtus_set3_sub0.6_run3.DAT mixtus_set3_sub0.6_run4.DAT mixtus_set3_sub0.6_run5.DAT mixtus_set4_sub0.6_run1.DAT mixtus_set4_sub0.6_run2.DAT mixtus_set4_sub0.6_run3.DAT mixtus_set4_sub0.6_run4.DAT mixtus_set4_sub0.6_run5.DAT mixtus_set5_sub0.6_run1.DAT mixtus_set5_sub0.6_run2.DAT mixtus_set5_sub0.6_run3.DAT mixtus_set5_sub0.6_run4.DAT mixtus_set5_sub0.6_run5.DAT impatiens_set1_sub0.6_run1.DAT impatiens_set1_sub0.6_run2.DAT impatiens_set1_sub0.6_run3.DAT impatiens_set1_sub0.6_run4.DAT impatiens_set1_sub0.6_run5.DAT impatiens_set2_sub0.6_run1.DAT impatiens_set2_sub0.6_run2.DAT impatiens_set2_sub0.6_run3.DAT impatiens_set2_sub0.6_run4.DAT impatiens_set2_sub0.6_run5.DAT impatiens_set3_sub0.6_run1.DAT impatiens_set3_sub0.6_run2.DAT impatiens_set3_sub0.6_run3.DAT impatiens_set3_sub0.6_run4.DAT impatiens_set3_sub0.6_run5.DAT impatiens_set4_sub0.6_run1.DAT impatiens_set4_sub0.6_run2.DAT impatiens_set4_sub0.6_run3.DAT impatiens_set4_sub0.6_run4.DAT impatiens_set4_sub0.6_run5.DAT impatiens_set5_sub0.6_run1.DAT impatiens_set5_sub0.6_run2.DAT impatiens_set5_sub0.6_run3.DAT impatiens_set5_sub0.6_run4.DAT impatiens_set5_sub0.6_run5.DAT 
)

# Run COLONY on the file for this array task
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}
./colony2s.gnu.out IFN:$INPUT
