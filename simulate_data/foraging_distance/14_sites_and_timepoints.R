######################################################################
# Simulate data with multiple landscapes and timepoints
######################################################################
# Started August 13, 2025
# With help from MB!


######################################################################
# Setup
######################################################################
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging")

par(family="serif", las=1, bty="l",
    cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 1))

library(rstan)
rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")


util <- new.env()
source('simulate_data/from_mike/mcmc_analysis_tools_rstan.R', local=util)
source('simulate_data/from_mike/mcmc_visualization_tools.R', local=util)

######################################################################
# Simulate And Explore Data
######################################################################

simu <- stan(file="simulate_data/src/sq_simu.stan",
             algorithm="Fixed_param", 
             seed=8438338, warmup=0, iter=1, chains=1, refresh=1)
delta <- rstan::extract(simu)$delta[1,,]
rho <- rstan::extract(simu)$rho[1]
trap_pos <- rstan::extract(simu)$trap_pos[1,,]
yobs <- rstan::extract(simu)$yobs[1,,]
sq_dist <- rstan::extract(simu)$sq_dist[1,,]

# Visualize detector array and source locations
par(mfrow=c(1, 1))

plot(delta[,1], delta[,2],
     col=util$c_light, pch=16, cex=0.5,
     xlab="x", xlim=c(0, 4500), 
     ylab="y", ylim=c(0, 3000))

points(trap_pos[,1], trap_pos[,2],
       col=util$c_dark, pch=16, cex=0.75)


par(mfrow=c(3, 3), mar = c(1, 1, 1, 1))

for (s in 1:9) {
  plot(trap_pos[,1], trap_pos[,2],
       col="black", pch=16, cex=sq_dist[s,] / 10^6,
       xlim=c(0, 4500), ylim=c(0, 3000))
  
  points(delta[s,1], delta[s,2],
         col=util$c_mid, pch=16, cex=0.8)
}



par(mfrow=c(10, 10), mar = c(1, 1, 1, 1))

for (s in 1:100) {
  plot(trap_pos[,1], trap_pos[,2],
       col="black", pch=16, cex=yobs[s,] / 2,
       xlim=c(0,4500), ylim=c(0,3000))
  
  points(delta[s,1], delta[s,2],
         col=util$c_mid, pch=16, cex=0.8)

}

tbl <- table(rowSums(yobs)[rowSums(yobs) > 0])
prop_tbl <- prop.table(tbl)  # gives proportions
barplot(prop_tbl, ylab = "Proportion", xlab = "Count category")

