# Compile summary results
# (e.g., output file didn't save but don't want to rerun all of the simulations)


param_grid = readRDS("simulate_data/methods_comparison/param_grid.rds")



for (task in 1:400){
  # load in data
  inputfilepath <- sprintf("simulate_data/methods_comparison/data/sim_result_%03d", task)
  yobs = readRDS(paste(inputfilepath, "/yobs.RDS", sep = ""))
  colony_data = readRDS(paste(inputfilepath, "/colonydata.RDS", sep = ""))
  trap_data = readRDS(paste(inputfilepath, "/trapdata.RDS", sep = ""))
  
  # add values to save
  nonzero = yobs[rowSums(yobs) > 0,]
  zero = yobs[rowSums(yobs) ==0,]
  param_grid$counts[param_grid$task_id == task] = list(rowSums(nonzero))
  param_grid$num_unobserved[param_grid$task_id == task] = nrow(zero)
  param_grid$true_average_foraging[param_grid$task_id == task] = mean(colony_data$foraging_range)
  param_grid$true_sd_foraging[param_grid$task_id == task] = sd(colony_data$foraging_range)
  
}

saveRDS(param_grid, "simulate_data/methods_comparison/output.rds")
