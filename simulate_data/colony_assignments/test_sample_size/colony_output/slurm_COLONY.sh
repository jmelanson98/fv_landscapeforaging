#!/bin/bash
#SBATCH --job-name=colony_run
#SBATCH --array=0-119
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=18:00:00
#SBATCH --output=logs/colony_%A_%a.out
#SBATCH --error=logs/colony_%A_%a.err

# Load any needed modules
module load gcc

# Your list of _no_exclusion.DAT files (replace with real names or glob dynamically)
FILES=(mixtus_set3_sub1_exclusion.DAT mixtus_set3_sub0.8_exclusion.DAT mixtus_set3_sub0.6_exclusion.DAT mixtus_set3_sub0.4_exclusion.DAT mixtus_set3_sub0.2_exclusion.DAT mixtus_set3_sub0.1_exclusion.DAT 
mixtus_set2_sub1_exclusion.DAT mixtus_set2_sub0.8_exclusion.DAT mixtus_set2_sub0.6_exclusion.DAT mixtus_set2_sub0.4_exclusion.DAT mixtus_set2_sub0.2_exclusion.DAT mixtus_set2_sub0.1_exclusion.DAT 
mixtus_set1_sub1_exclusion.DAT mixtus_set1_sub0.8_exclusion.DAT mixtus_set1_sub0.6_exclusion.DAT mixtus_set1_sub0.4_exclusion.DAT mixtus_set1_sub0.2_exclusion.DAT mixtus_set1_sub0.1_exclusion.DAT 
mixtus_set4_sub1_exclusion.DAT mixtus_set4_sub0.8_exclusion.DAT mixtus_set4_sub0.6_exclusion.DAT mixtus_set4_sub0.4_exclusion.DAT mixtus_set4_sub0.2_exclusion.DAT mixtus_set4_sub0.1_exclusion.DAT 
mixtus_set5_sub1_exclusion.DAT mixtus_set5_sub0.8_exclusion.DAT mixtus_set5_sub0.6_exclusion.DAT mixtus_set5_sub0.4_exclusion.DAT mixtus_set5_sub0.2_exclusion.DAT mixtus_set5_sub0.1_exclusion.DAT 
impatiens_set3_sub1_exclusion.DAT impatiens_set3_sub0.8_exclusion.DAT impatiens_set3_sub0.6_exclusion.DAT impatiens_set3_sub0.4_exclusion.DAT impatiens_set3_sub0.2_exclusion.DAT impatiens_set3_sub0.1_exclusion.DAT 
impatiens_set2_sub1_exclusion.DAT impatiens_set2_sub0.8_exclusion.DAT impatiens_set2_sub0.6_exclusion.DAT impatiens_set2_sub0.4_exclusion.DAT impatiens_set2_sub0.2_exclusion.DAT impatiens_set2_sub0.1_exclusion.DAT 
impatiens_set1_sub1_exclusion.DAT impatiens_set1_sub0.8_exclusion.DAT impatiens_set1_sub0.6_exclusion.DAT impatiens_set1_sub0.4_exclusion.DAT impatiens_set1_sub0.2_exclusion.DAT impatiens_set1_sub0.1_exclusion.DAT 
impatiens_set1_sub4_exclusion.DAT impatiens_set4_sub0.8_exclusion.DAT impatiens_set4_sub0.6_exclusion.DAT impatiens_set4_sub0.4_exclusion.DAT impatiens_set4_sub0.2_exclusion.DAT impatiens_set4_sub0.1_exclusion.DAT 
impatiens_set1_sub5_exclusion.DAT impatiens_set5_sub0.8_exclusion.DAT impatiens_set5_sub0.6_exclusion.DAT impatiens_set5_sub0.4_exclusion.DAT impatiens_set5_sub0.2_exclusion.DAT impatiens_set5_sub0.1_exclusion.DAT 
mixtus_set3_sub1_no_exclusion.DAT mixtus_set3_sub0.8_no_exclusion.DAT mixtus_set3_sub0.6_no_exclusion.DAT mixtus_set3_sub0.4_no_exclusion.DAT mixtus_set3_sub0.2_no_exclusion.DAT mixtus_set3_sub0.1_no_exclusion.DAT 
mixtus_set2_sub1_no_exclusion.DAT mixtus_set2_sub0.8_no_exclusion.DAT mixtus_set2_sub0.6_no_exclusion.DAT mixtus_set2_sub0.4_no_exclusion.DAT mixtus_set2_sub0.2_no_exclusion.DAT mixtus_set2_sub0.1_no_exclusion.DAT 
mixtus_set1_sub1_no_exclusion.DAT mixtus_set1_sub0.8_no_exclusion.DAT mixtus_set1_sub0.6_no_exclusion.DAT mixtus_set1_sub0.4_no_exclusion.DAT mixtus_set1_sub0.2_no_exclusion.DAT mixtus_set1_sub0.1_no_exclusion.DAT 
mixtus_set4_sub1_no_exclusion.DAT mixtus_set4_sub0.8_no_exclusion.DAT mixtus_set4_sub0.6_no_exclusion.DAT mixtus_set4_sub0.4_no_exclusion.DAT mixtus_set4_sub0.2_no_exclusion.DAT mixtus_set4_sub0.1_no_exclusion.DAT 
mixtus_set5_sub1_no_exclusion.DAT mixtus_set5_sub0.8_no_exclusion.DAT mixtus_set5_sub0.6_no_exclusion.DAT mixtus_set5_sub0.4_no_exclusion.DAT mixtus_set5_sub0.2_no_exclusion.DAT mixtus_set5_sub0.1_no_exclusion.DAT 
impatiens_set3_sub1_no_exclusion.DAT impatiens_set3_sub0.8_no_exclusion.DAT impatiens_set3_sub0.6_no_exclusion.DAT impatiens_set3_sub0.4_no_exclusion.DAT impatiens_set3_sub0.2_no_exclusion.DAT impatiens_set3_sub0.1_no_exclusion.DAT 
impatiens_set2_sub1_no_exclusion.DAT impatiens_set2_sub0.8_no_exclusion.DAT impatiens_set2_sub0.6_no_exclusion.DAT impatiens_set2_sub0.4_no_exclusion.DAT impatiens_set2_sub0.2_no_exclusion.DAT impatiens_set2_sub0.1_no_exclusion.DAT 
impatiens_set1_sub1_no_exclusion.DAT impatiens_set1_sub0.8_no_exclusion.DAT impatiens_set1_sub0.6_no_exclusion.DAT impatiens_set1_sub0.4_no_exclusion.DAT impatiens_set1_sub0.2_no_exclusion.DAT impatiens_set1_sub0.1_no_exclusion.DAT 
impatiens_set1_sub4_no_exclusion.DAT impatiens_set4_sub0.8_no_exclusion.DAT impatiens_set4_sub0.6_no_exclusion.DAT impatiens_set4_sub0.4_no_exclusion.DAT impatiens_set4_sub0.2_no_exclusion.DAT impatiens_set4_sub0.1_no_exclusion.DAT 
impatiens_set1_sub5_no_exclusion.DAT impatiens_set5_sub0.8_no_exclusion.DAT impatiens_set5_sub0.6_no_exclusion.DAT impatiens_set5_sub0.4_no_exclusion.DAT impatiens_set5_sub0.2_no_exclusion.DAT impatiens_set5_sub0.1_no_exclusion.DAT
)

# Run COLONY on the file for this array task
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}
./colony2s.gnu.out IFN:$INPUT
