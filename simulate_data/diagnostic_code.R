##### Some useful Stan diagnostic code #####
# Script Initiated: May 11, 2025
# Started by: Jenna Melanson
# This code is modified from various tidbits from Michael Betancourt (https://betanalpha.github.io)


# check distribution of divergent transitions in parameter space
partition <- partition_div(fit_cp)
div_params_cp <- partition[[1]]
nondiv_params_cp <- partition[[2]]

par(mar = c(4, 4, 0.5, 0.5))
plot(nondiv_params_cp$'theta[1]', log(nondiv_params_cp$tau),
     col=c_dark, pch=16, cex=0.8, xlab="theta.1", ylab="log(tau)",
     xlim=c(-20, 50), ylim=c(-6,4))
points(div_params_cp$'theta[1]', log(div_params_cp$tau),
       col="green", pch=16, cex=0.8)


stan_rhat
stan_ess
traceplot(exfit, pars = "sigma", inc_warmup = TRUE)
