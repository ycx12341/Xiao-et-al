### Gradient matching scheme ########
### Author: Yunchen Xiao & Len Thomas ###########

#Source companion functions
source("PDE_GradientMatching_Functions.r")

#Load the package "tictoc" in order to measure the computational time. 
library(tictoc) 
#Load "readr" package for writing results - just check it's loaded
#(it's actually used first inside the parallel routine)
library(readr)
#Load packages for running the simulation in parallel
library(doParallel) 
library(doRNG)


### Setup ####

#Define simulation parameters
#Measurement error CV levels to run (if a scalar then just runs at this one level)
cv <- c(0.01, 0.025, 0.05, 0.075, 0.10)
#Measurement error distribution
dist = "gamma"
#Number of simulations at each level of CV
n.sims <- 200
#Whether to save simulation outputs to file or not
# (useful as they are time-consuming to run, and this
# allows further post-processing on results)
save.sims <- TRUE
save.sims.dir <- "SimRes"
#Set random number seed so results are reproducible
# (random number set within each level of CV to make any repeat runs of a 
#  single CV level easier)
rn.seed <- 874513
#Number of parallel threads to run on
n.threads <- 25 #detectCores() - 1

# Define model parameters
dn <- 0.01
gamma <- 0.05
eta <- 10
dm <- 0.01
alpha <- 0.1
rn <- 5
# This parameter not included in the optimization
beta <- 0
# Make a vector to store the true values
true.values <- c(dn, gamma, rn, eta, dm, alpha)
names(true.values) <- c("dn", "gamma", "rn", "eta", "dm", "alpha")

# Define 1D dimensionless space points
n.x11 <- 80
max.x11 <- 1
x11 <- seq(0, max.x11, length = n.x11)

# Define time discretization and max time
dt <- 0.001
max.t <- 10

# Set initial conditions
eps <- 0.01
n0 <- rep(0, n.x11)
for (i in 1:n.x11) {
  if (x11[i] <= 0.25) {
    n0[i] <- exp(-x11[i] ^ 2 / eps)
  } else {
    n0[i] <- 0
  }
}

# Initial values
n <- n0
f0 <- 1-0.5*n0
f <- f0
m0 <- 0.5*n0
m <- m0

#Generate reference dataset
ref.data <- generate.reference.data(n.x11, max.x11, dt, max.t,
                                    dn, gamma, eta, dm, alpha, rn, beta, n0, f0, m0)
data.dim.row <- length(ref.data$n[, 1])
data.dim.col <- length(ref.data$n[1, ])
#Truncate the reference data set to avoid generating gamma errors at data values of 0 
ref.data.trun <- list(n = ref.data$n[2:(data.dim.row - 1), 2:(data.dim.col - 1)], 
                      f = ref.data$f[2:(data.dim.row - 1), 2:(data.dim.col - 1)], 
                      m = ref.data$m[2:(data.dim.row - 1), 2:(data.dim.col - 1)])

#Use start values from manuscript
start.values <- c(1E-2, 1.28E-1, 6.30, 1.30E1, 1.70E-2, 1.30E-1)

#Record degrees of freedom for sd calculations
n.data.trun <- sum(apply(sapply(ref.data.trun, dim), 2, prod))
n.params <- length(start.values)
df <- n.data.trun - n.params


### Run simulations ####

if(save.sims) {
  if(!dir.exists(save.sims.dir)) dir.create(save.sims.dir)
}

cl <- makeCluster(n.threads) 
registerDoParallel(cl)
#Set parallel seed - once = TRUE means it's set each time foreach is called
registerDoRNG(rn.seed, once = TRUE)
#Start the timer
tic()
for(i in 1:length(cv)) {
  #Run the simulation in parallel
  ests <- foreach (sim = 1:n.sims, .combine = rbind) %dopar% {
    #Introduce random error
    pert.data <- perturb.reference.data(ref.data.trun, cv = cv[i], distribution = dist)
    
    #Obtain gradient approximations
    grads <- approximate.gradients(pert.data, x11, max.t, distribution = dist)
    
    #Estimate parameter values and associated sds
    res <- optim(start.values, calculate.sse, grads = grads, hessian = TRUE)
    par.ests <- res$par
    sd.ests <- estimate.sd(res, df)
    
    #Save sim results to file, so they can readily be retrieved later
    if(save.sims)
      readr::write_rds(list(ref.data.trun = ref.data.trun, pert.data = pert.data, 
        grads = grads, cv = cv[i], dist = dist, res = res, par.ests = par.ests, 
        sd.ests = sd.ests), 
        path = paste0("./", save.sims.dir, "/cv", i, "_sim", sim, "_res.rds"))
  
    #Vector to return from the foreach
    c(par.ests, sd.ests)
  }
  par.ests <- ests[, 1:length(true.values)]
  sd.ests <- ests[, 1:length(true.values) + length(true.values)]
  colnames(par.ests) <- colnames(sd.ests) <- names(true.values)
  if(save.sims) 
    write_rds(list(true.values = true.values, par.ests = par.ests, sd.ests = sd.ests),
              paste0("./", save.sims.dir, "/cv", i, "_sim_res.rds"))
}
#Stop the timer
toc()
#Stop the cluster
stopCluster(cl)


