###### Fit foraging distance model to real data!
### November 14, 2025
### J. Melanson

##### Load packages #####
library(rstan)
library(matrixStats)
library(sp)
library(gstat)
library(ggplot2)
library(reshape2)
library(raster)
library(rasterVis)
library(parallel)
library(future)
library(furrr)
library(dplyr)
library(tidyr)
library(gridExtra)
library(tibble)
library(sf)

##### Set Environment #####
#setwd("/Users/jenna1/fv_landscapeforaging")
#bombus_path = "/Users/jenna1/Documents/UBC/bombus_project"

setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "/home/melanson/projects/def-ckremen/melanson"
source("src/analysis_functions.R")


# Load in data
mixtus_sibs2022 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mixtus_sibs2023 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
impatiens_sibs2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
impatiens_sibs2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")
specs2022 = read.csv(paste0(bombus_path, "/raw_data/2022specimendata.csv"))
specs2023 = read.csv(paste0(bombus_path, "/raw_data/2023specimendata.csv"))
effort2022 = read.csv(paste0(bombus_path, "/raw_data/2022sampledata.csv"))
effort2023 = read.csv(paste0(bombus_path, "/raw_data/2023sampledata.csv"))
samplepoints = read.csv(paste0(bombus_path, "/raw_data/allsamplepoints.csv"), header = FALSE)



# Prepare data for Stan (with function)
mix22_stan = prep_stan_simpleforaging(mixtus_sibs2022,
                                      specs2022,
                                      effort2022,
                                      samplepoints)
mix23_stan = prep_stan_simpleforaging(mixtus_sibs2023,
                                      specs2023,
                                      effort2023,
                                      samplepoints)
imp22_stan = prep_stan_simpleforaging(impatiens_sibs2022,
                                      specs2022,
                                      effort2022,
                                      samplepoints)
imp23_stan = prep_stan_simpleforaging(impatiens_sibs2023,
                                      specs2023,
                                      effort2023,
                                      samplepoints)



#select stan model to fit
stanfile = "models/simple_multinomial.stan"


# Get task ID from slurm manager
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
data_list = list(mix22_stan, mix23_stan, imp22_stan, imp23_stan)
data = data_list[[task_id]]

#fit and save model
stanFit = stan(file = stanfile,
               data = data, seed = 5838299,
               chains = 4, cores = 4,
               control = list(max_treedepth = 15),
               iter = 4000,
               verbose = TRUE)
saveRDS(stanFit, paste0("analysis/foragingmodel_", task_id, ".rds"))


#Plot the posteriors of some colonies
plot_list = list()
numplots = 3
legends = list()

for (i in 1:numplots){
  delta_draws = cbind(x = rstan::extract(stanFitED, pars = "delta_x")$delta[, i],
                      y = rstan::extract(stanFitED, pars = "delta_y")$delta[, i])
  traps_temp = full_join(traps_m_2023, filled_counts[filled_counts$sibshipID ==i,])
  traps_temp$count[is.na(traps_temp$count)] = 0
  traps_temp$trap_x = traps_temp$trap_x/1000
  traps_temp$trap_y = traps_temp$trap_y/1000
  
  p = ggplot(delta_draws, aes(x = x, y = y)) +
    geom_density_2d_filled(alpha = 0.8) +
    
    #plot trap locations / sizes / quality
    geom_point(data = traps_temp, aes(x = trap_x, y = trap_y, size = count, colour = "red")) +
    scale_size_continuous(limits = c(0,10), range = c(1, 5)) +
    
    #miscellaneous
    labs(title = paste("Colony", i), 
         size = "Number of Captures",
         level = "Colony Posterior") +
    guides(colour = "none") +
    #xlim(c(1220, 1225)) +
    #ylim(c(455,460)) +
    coord_equal() +
    theme_bw()
  
  # save legend
  g <- ggplotGrob(p)
  legend_index <- which(g$layout$name == "guide-box-right")
  legend <- g$grobs[[legend_index]]
  
  # remove legend from plot
  p <- p + theme(legend.position = "none")
  
  #save plot
  plot_list[[i]] = p
  legends[[1]] = legend
}

fig = grid.arrange(grobs = plot_list, ncol = 3)
fig = grid.arrange(fig, legends[[1]], ncol = 2, widths = c(4,1))
ggsave(paste0("analysis/foragingmodel_", task_id, ".jpg"), fig, height = 3000, width = 4000, units = "px")

print("Done!")
