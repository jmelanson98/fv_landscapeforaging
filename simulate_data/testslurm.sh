#!/bin/bash
##SBATCH --job-name=test
#SBATCH --output=logs/test_%A_%a.out
#SBATCH --error=logs/test_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=00:10:00

module load StdEnv/2023 gcc r/4.3.1 gdal proj

# Define your personal R library path
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.3

#!/bin/bash

Rscript - <<EOF
library(cmdstanr)
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
set_cmdstan_path("/home/melanson/projects/def-ckremen/melanson/cmdstan")

mod <- cmdstan_model(
  "/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging/models/test.stan",
  force_recompile = TRUE
)

saveRDS(mod, "compiled_test_model.rds")
EOF


