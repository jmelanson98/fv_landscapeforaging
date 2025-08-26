#!/bin/bash
#SBATCH --job-name=colony_run
#SBATCH --array=0-23
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/colony_%A_%a.out
#SBATCH --error=logs/colony_%A_%a.err

# Load any needed modules
module load gcc

# Your list of .DAT files (replace with real names or glob dynamically)
FILES=(mixtus_pp0_poly.DAT mixtus_pp0.DAT
mixtus_pp0.2_poly.DAT mixtus_pp0.2.DAT
mixtus_pp0.4_poly.DAT mixtus_pp0.4.DAT
mixtus_pp0.6_poly.DAT mixtus_pp0.6.DAT
mixtus_pp0.8_poly.DAT mixtus_pp0.8.DAT
mixtus_pp1_poly.DAT mixtus_pp1.DAT impatiens_pp0_poly.DAT impatiens_pp0.DAT
impatiens_pp0.2_poly.DAT impatiens_pp0.2.DAT
impatiens_pp0.4_poly.DAT impatiens_pp0.4.DAT
impatiens_pp0.6_poly.DAT impatiens_pp0.6.DAT
impatiens_pp0.8_poly.DAT impatiens_pp0.8.DAT
impatiens_pp1_poly.DAT impatiens_pp1.DAT
)

# Run COLONY on the file for this array task
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}
./colony2s.gnu.out IFN:$INPUT
