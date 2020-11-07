### Calculation of Bootstrap SD from Gradient matching scheme########
### Author: Yunchen Xiao & Len Thomas ###########

#Calculates bootstrap SD from graident matching -- run
# PDE_GradientMatching_Main.r first to generate results
# that are then bootstrapped with this code.

#Source companion functions
source("PDE_GradientMatching_Functions.r")

#readr used to read and write results to file
library(readr)
#Load packages for running the simulation in parallel
library(doParallel) 
library(doRNG)

### Setup ###

#Directory to read sim results from and write bootstrap results to
save.sims.dir <- "SimRes1"
#Measurement error CVs that were run
cv <- c(0.01, 0.025, 0.05, 0.075, 0.10)
n.cvs <- length(cv)

#Number of bootstrap replicates to calculate SD
B <- 200

#Start values to use for parameters in optim
start.values <- c(1E-2, 1.28E-1, 6.30, 1.30E1, 1.70E-2, 1.30E-1)
n.pars <- length(start.values)

#Set random number seed so results are reproducible
# (random number set within each level of CV to make any repeat runs of a 
#  single CV level easier)
rn.seed <- 762905
#Number of parallel threads to run on
n.threads <- 25 #detectCores() - 1


### Bootstrap ###

dim.grads <- dim(read_rds(paste0(save.sims.dir, "/cv1_sim1_res.rds"))$grads[[1]])
n.data <- prod(dim.grads)
n.sims <- nrow(read_rds(paste0(save.sims.dir, "/cv1_sim_res.rds"))$par.ests)
boot.res <- matrix(NA, B, n.pars)

cl <- makeCluster(n.threads) 
registerDoParallel(cl)
#Set parallel seed - once = TRUE means it's set each time foreach is called
registerDoRNG(rn.seed, once = TRUE)
for(i in 1:n.cvs) {
  ests <- foreach (k = 1:n.sims, .combine = rbind) %dopar% {
    sim <- readr::read_rds(paste0(save.sims.dir, "/cv", i, "_sim", k, "_res.rds"))
    boot.grads <- sim$grads
    for (b in 1:B) {
      new.dat.ind <- sample.int(n.data, replace = TRUE)
      for(l in 1:length(sim$grads))
        boot.grads[[l]] <- matrix(as.numeric(sim$grads[[l]])[new.dat.ind], 
          dim.grads[1], dim.grads[2])
      boot.res[b, ] <- optim(start.values, calculate.sse, grads = boot.grads)$par
    }
    readr::write_rds(boot.res, paste0(save.sims.dir, "/cv", i, "_sim", k, "_boot_res.rds"))
    boot.mean <- apply(boot.res, 2, mean)
    boot.sd <- apply(boot.res, 2, sd)
    #Vector to return from the foreach
    c(boot.mean, boot.sd)
  }
  boot.mean <- ests[, 1:n.pars]
  boot.sd <- ests[, 1:n.pars + n.pars]
  readr::write_rds(list(boot.mean = boot.mean, boot.sd = boot.sd), 
            paste0(save.sims.dir, "/cv", i, "_boot_res.rds"))
}
#Stop the cluster
stopCluster(cl)



