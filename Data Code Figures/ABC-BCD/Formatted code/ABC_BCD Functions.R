#### ABC-BCD scheme - Functions ########
### AuthorS: Yunchen Xiao & Len Thomas ###########

generate.init.paras <- function(x, rng.seed) {
# Purpose: Generates the initial parameters set.
# Inputs: 
# x:  Number of parameter vectors being generated.
# Output:
# init.paras: x rows of initial parameter vectors, combined
# into a matrix.
# Implementation note: 
# set.seed() is used in the function to generate reproducible
# initial parameter set.
  
  # Set seed for the random number generator.
  if (!is.null(rng.seed)){
    set.seed(rng.seed)
  }
  
  # Sample parameter values from their corresponding initial distribution.
  dn <- runif(x, 0.000069, 0.02)
  gamma <- runif(x, 0.005, 0.26)
  eta <- runif(x, 7, 18)
  dm <- runif(x, 0.0001, 0.033)
  alpha <- runif(x, 0.07, 0.18)
  rn <- runif(x, 3.5, 9)
  
  # Form the initial parameter set.
  init.paras <- cbind(dn, gamma, eta, dm, alpha, rn)
  return(init.paras)  
}

generate.ref.ss <- function(ref.paras) {
# Purpose: Generates summary statistics based on the reference dataset.
# Inputs: 
# ref.paras: Reference parameter values
# Output:
# mean.var: Mean and variance of each time series in the reference dataset.
  
  # Space
  l1 <- 0
  l2 <- 1
  x11 <- seq(l1, l2, length = 80)
  n.x11 <- length(x11)
  h <- x11[2] - x11[1]
  
  # Time
  T <- 10
  dt <- 0.001
  t <- seq(0, T, by = 0.001)
  
  # Reference parameter values
  dn <- ref.paras[1]
  gamma <- ref.paras[2]
  ita <- ref.paras[4]
  dm <- ref.paras[5]
  alpha <- ref.paras[6]
  r <- ref.paras[3]
  
  # Parameter values excluded from the inference
  beta <- 0
  eps <- 0.01
  
  # Initial condition
  n0 <- rep(0, length(x11))
  for (i in 1:length(x11)) {
    if (x11[i] <= 0.25) {
      n0[i] <- exp(-x11[i] ^ 2 / eps)
    } else {
      n0[i] <- 0
    }
  }
  f0 <- 1 - 0.5*n0
  m0 <- 0.5 * n0
  
  # Initial densities for the numerical solver
  n <- n0
  f <- f0
  m <- m0
    
  p <- 1
    
  # Initialise an empty vector to store the reference dataset.
  # n - tumour cells
  # f - ECM
  # m - MDE

  res <- matrix(0, nrow = 3 * T, ncol = n.x11)
    
  while(p * dt <= T) {
    f[2:(n.x11 - 1)] <- -ita * dt * m[2:(n.x11 - 1)] * f[2:(n.x11 - 1)] + 
      f[2:(n.x11 - 1)]
      
    m[2:(n.x11 - 1)] <- dm * (m[1:(n.x11 - 2)] + m[3:n.x11] - 
      2 * m[2:(n.x11 - 1)]) * dt / (h ^ 2) +alpha * n[2:(n.x11 - 1)] * dt -  
      beta * m[2:(n.x11 - 1)] * dt + m[2:(n.x11 - 1)]
      
    n[2:(n.x11 - 1)] <- dn * (n[1:(n.x11 - 2)] + n[3:n.x11] - 2 * 
      n[2:(n.x11 - 1)]) * dt / (h ^ 2) -  gamma * (n[3:n.x11] - 
      n[2:(n.x11 - 1)]) * (f[3:n.x11]-f[2:(n.x11 - 1)]) * dt / (h ^ 2) - 
      gamma * n[2:(n.x11 - 1)] * (f[1:(n.x11 - 2)] + f[3:n.x11] - 
      2 * f[2:(n.x11 - 1)]) * dt / (h ^ 2) + 
      r * (1 - f[2:(n.x11 - 1)] - n[2:(n.x11 - 1)]) * n[2:(n.x11 - 1)]* dt + 
      n[2:(n.x11 - 1)]
      
    # No flux boundary condition
    n[1] <- n[2]
    n[n.x11] <- n[n.x11 - 1]
      
    f[1] <- f[2]
    f[n.x11] <- f[n.x11 - 1]
      
    m[1] <- m[2]
    m[n.x11] <- m[n.x11 - 1]
      
    # Save the results at integer time steps
    if(p %% (1 / dt) == 0) {
      ind.cohort <- p / (1 / dt)
      res[(ind.cohort * 3 - 2), ] <- n
      res[(ind.cohort * 3 - 1), ] <- f
      res[ind.cohort * 3, ] <- m
    }
      
    p <- p + 1
  }
  
  # Rearrange the reference dataset in the order of:
  # tumour cells - ECM - MDE
  res.arr <- matrix(0,nrow = 3 * T, ncol = n.x11)
    
  tc.ind <- seq(1, (3 * T - 2), by = 3)
  ecm.ind <- seq(2, (3 * T - 1), by = 3)
  mde.ind <- seq(3, (3 * T), by = 3)
    
  res.arr[1:T, ] = res[tc.ind, ]
  res.arr[(T + 1):(2 * T), ] = res[ecm.ind, ]
  res.arr[(2 * T + 1):(3 * T), ] = res[mde.ind, ]
  
  # Set an empty matrix to store the summary statistics.  
  mean.var = matrix(0, nrow = n.x11 * 3, ncol = 2) 
    
  for (i in 1:n.x11) {
    mean.var[i, 1] <- mean(res.arr[1:T, i])
    mean.var[i, 2] <- var(res.arr[1:T, i])
      
    mean.var[i + n.x11, 1] <- mean(res.arr[(T + 1):(2 * T), i])
    mean.var[i + n.x11, 2] <- var(res.arr[(T + 1):(2 * T), i])
    mean.var[i + (2 * n.x11), 1] <- mean(res.arr[(2 * T + 1):(3 * T), i])
    mean.var[i + (2 * n.x11), 2] <- var(res.arr[(2 * T + 1):(3 * T), i])
  }
  
  return(mean.var)
}


bcd <- function(ref.paras, paras, ind, print.bcd) {
# Purpose: Calculates each parameter vector's Bhattacharyya distance to the 
# reference parameter vector.
# Inputs: 
# ref.paras: Reference parameter values.
# paras: Parameter sets.
# ind: An indicator suggests which trial the input parameter set is related to.
# "ecm" - the first trial, ECM-related rounds.
# "ecm_mde" - the second trial, ECM+MDE_related rounds.
# "all_3" - the final trial, the Bhattacharyya distance across all 3 equations 
# are calculated. 
# Output:
# bcd.tab: Each parameter vector's Bhattacharyya distance to the reference
# parameter vector.

  # Reference summary statistics
  mean.var.obs <- generate.ref.ss(ref.paras)
  
  # Space
  l1 <- 0
  l2 <- 1
  x11 <- seq(l1, l2, length = 80)
  n.x11 <- length(x11)
  h <- x11[2] - x11[1]
  
  # Time
  T <- 10
  dt <- 0.001
  t <- seq(0, T, by = 0.001)
  
  # Parameter values
  dn <- paras[, 1]
  gamma <- paras[, 2]
  ita <- paras[, 3]
  dm <- paras[, 4]
  alpha <- paras[, 5]
  r <- paras[, 6]
  
  # Parameter values excluded from the inference
  beta <- 0
  eps <- 0.01
  
  # Initial condition
  n0 <- rep(0, length(x11))
  for (i in 1:length(x11)) {
    if (x11[i] <= 0.25) {
      n0[i] <- exp(-x11[i] ^ 2 / eps)
    } else {
      n0[i] <- 0
    }
  }
  f0 <- 1-0.5*n0
  m0 <- 0.5*n0
  
  # Initialise an empty vector, used later to store the 
  # Bhattacharyya distances.
  
  bcd.tab <- matrix(0, nrow = length(paras[, 1]), ncol = 2)
  
  for (q in 1:length(paras[,1])) {
    
    # Initial densities for the numerical solver.
    
    n <- n0
    f <- f0
    m <- m0
    
    p <- 1

    res <- matrix(0, nrow = T * 3, ncol = n.x11)
    
    while(p * dt <= T) {
      f[2:(n.x11 - 1)] <- -ita[q] * dt * m[2:(n.x11 - 1)] * f[2:(n.x11 - 1)] + 
        f[2:(n.x11 - 1)]
      
      m[2:(n.x11 - 1)] <- dm[q] * (m[1:(n.x11 - 2)] + m[3:n.x11] - 
        2 * m[2:(n.x11 - 1)]) * dt / (h ^ 2) +alpha[q] * n[2:(n.x11 - 1)] * dt -  
        beta * m[2:(n.x11 - 1)] * dt + m[2:(n.x11 - 1)]
      
      n[2:(n.x11 - 1)] <- dn[q] * (n[1:(n.x11 - 2)] + n[3:n.x11] - 2 * 
        n[2:(n.x11 - 1)]) * dt / (h ^ 2) -  gamma[q] * (n[3:n.x11] - 
        n[2:(n.x11 - 1)]) * (f[3:n.x11]-f[2:(n.x11 - 1)]) * dt / (h ^ 2) - 
        gamma[q] * n[2:(n.x11 - 1)] * (f[1:(n.x11 - 2)] + f[3:n.x11] - 
        2 * f[2:(n.x11 - 1)]) * dt / (h ^ 2) + r[q] * (1 - f[2:(n.x11 - 1)] - 
        n[2:(n.x11 - 1)]) * n[2:(n.x11 - 1)] * dt + n[2:(n.x11 - 1)]
      
      # No flux boundary condition
      n[1] <- n[2]
      n[n.x11] <- n[n.x11 - 1]
      
      f[1] <- f[2]
      f[n.x11] <- f[n.x11 - 1]
      
      m[1] <- m[2]
      m[n.x11] <- m[n.x11 - 1]
      
      # Save the results at each integer time steps
      if(p %% (1/dt) == 0) {
        ind.cohort <- p / (1 / dt)
        res[(ind.cohort * 3 - 2), ] <- n
        res[(ind.cohort * 3 - 1), ] <- f
        res[ind.cohort * 3, ] <- m
      }
      
      p <- p + 1
    }
    
    # Rearrange the summary statistics in the same order as the reference
    # summary statistics.
    
    res.arr <- matrix(0,nrow = 3 * T, ncol = n.x11)
    
    tc.ind <- seq(1, (3 * T - 2), by = 3)
    ecm.ind <- seq(2, (3 * T - 1), by = 3)
    mde.ind <- seq(3, (3 * T), by = 3)
    
    res.arr[1:T, ] = res[tc.ind,]
    res.arr[(T + 1):(2 * T), ] = res[ecm.ind,]
    res.arr[(2 * T + 1):(3 * T), ] = res[mde.ind,]
    
    mean.var = matrix(0, nrow = length(mean.var.obs[, 1]), 
     ncol = length(mean.var.obs[1, ])) 
    
    for (i in 1:n.x11) {
      mean.var[i, 1] <- mean(res.arr[1:T, i])
      mean.var[i, 2] <- var(res.arr[1:T, i])
      
      mean.var[i + n.x11, 1] <- mean(res.arr[(T + 1):(2 * T), i])
      mean.var[i + n.x11, 2] <- var(res.arr[(T + 1):(2 * T), i])
      
      mean.var[i + (2 * n.x11), 1] <- mean(res.arr[(2 * T + 1):(3 * T),i])
      mean.var[i + (2 * n.x11), 2] <- var(res.arr[(2 * T + 1):(3 * T),i])
    }
    
    # Initialise an empty vector to store the Bhattacharyya distance of each time
    # series in the simulated dataset.
    # ind == "ecm": store the results from the second time series compartment.
    # ind == "ecm_mde": store the results from second and the third time 
    # series compartment.
    # ind == "all_3": store the results from all time series.
    bcd.vec <- vector()
    
    if (ind == "ecm") {
      for (j in (n.x11 + 1):(2 * n.x11)) {
        bcd <- 0.25 * log(0.25 * ((mean.var[j, 2] / mean.var.obs[j, 2]) + 
              (mean.var.obs[j, 2] / mean.var[j, 2]) + 2)) + 0.25 * 
              (((mean.var[j, 1] - mean.var.obs[j, 1]) ^ 2) / (mean.var.obs[j, 2] + 
              mean.var[j, 2]))
        bcd.vec <- c(bcd.vec,bcd)
      } 
    } else if (ind == "ecm_mde") {
      for (j in (n.x11 + 1):(3 * n.x11)) {
        bcd <- 0.25 * log(0.25 * ((mean.var[j, 2] / mean.var.obs[j, 2]) + 
              (mean.var.obs[j, 2] / mean.var[j, 2]) + 2)) + 0.25 * 
              (((mean.var[j, 1] - mean.var.obs[j, 1]) ^ 2) / (mean.var.obs[j, 2] + 
              mean.var[j, 2]))
        bcd.vec <- c(bcd.vec,bcd)
      } 
    } else if (ind == "all_3") {
      for (j in 1:(3 * n.x11)) {
        bcd <- 0.25 * log(0.25 * ((mean.var[j, 2] / mean.var.obs[j, 2]) + 
              (mean.var.obs[j, 2] / mean.var[j, 2]) + 2)) + 0.25 * 
              (((mean.var[j, 1] - mean.var.obs[j, 1]) ^ 2) / (mean.var.obs[j, 2] + 
              mean.var[j, 2]))
        bcd.vec <- c(bcd.vec,bcd)
      } 
    }
    
    # Locate the invalid terms and set them to be 0
    inv.index <- which(bcd.vec == "Inf")
    bcd.vec.2 <- bcd.vec
    bcd.vec.2[inv.index] <- 0
    
    # Sum up the Bhattacharyya distance across the time series in the corresponding
    # compartment.
    bcd.sum <- sum(bcd.vec.2)
    bcd.tab[q, ] <- c(q, bcd.sum)

    # Optional argument used to track the progress.
    if (print.bcd == TRUE) {
      print(c(q, bcd.sum))
    }
  }
  return(bcd.tab)
}

abc.bcd <- function(ss.mat,paras,rng.seed) {
# Purpose: Resample and perturb the parameter set of current round, 
# based on the Bhattacharyya distance results of each parameter vector.
# Inputs:
# ss.mat: Bhattacharyya distance results  
# paras: Parameter sets.
# Output:
# paras.nr.perturbed: Parameter set for the next round.
# Implementation note: 
# set.seed() is used in the function for reproducible results.

  if (!is.null(rng.seed)){
    set.seed(rng.seed)
  }

  # Set the Bhattacharya distance results to be a matrix 
  # instead of a data frame.
  ss.mat <- as.matrix(ss.mat) 
  
  # Locate the invalid terms in the Bhattacharyya distance results.
  invalid <- vector()
  
  for (j in 1:length(ss.mat[1, ])) {
    invalid.sep <- which(is.na(ss.mat[, j]) == TRUE)
    invalid <- c(invalid, invalid.sep)
  }  
  
  # Uniqueness checking, make sure each invalid term only appears once.
  invalid.index <- unique(invalid) 
  
  # Valid Bhattacharya distance results.
  if(length(invalid.index) == 0) {
    ss.mat.valid <- ss.mat
  } else {
    ss.mat.valid <- ss.mat[-invalid.index,]
  }
  
  # Weights of the valid Bhattacharyya distance results.
  wt <- 1/(ss.mat.valid[,length(ss.mat.valid[1,])]^(1/2)) 
  
  # Valid Bhattacharyya distance results + Weight
  ss.mat.valid.wt <- cbind(ss.mat.valid,wt)
  
  # Create an empty vector to store the resampling probabilities
  resamp.prob <- rep(0,length(ss.mat.valid.wt[,length(ss.mat.valid.wt[1,])])) 
  
  # Calculation of resampling probabilities, see the manuscript for more details.
  for (i in 1:length(ss.mat.valid.wt[, length(ss.mat.valid.wt[1, ])])) {
    if (ss.mat.valid.wt[i, length(ss.mat.valid.wt[1, ])] == 
        min(ss.mat.valid.wt[,length(ss.mat.valid.wt[1, ])])) {
      resamp.prob[i] <- 0
    } else if (ss.mat.valid.wt[i, length(ss.mat.valid.wt[1, ])] == 
               max(ss.mat.valid.wt[, length(ss.mat.valid.wt[1, ])])) {
      resamp.prob[i] <- 1
    } else {
      resamp.prob[i] <- (ss.mat.valid.wt[i, length(ss.mat.valid.wt[1, ])] - 
                         min(ss.mat.valid.wt[,length(ss.mat.valid.wt[1, ])])) / 
                         (max(ss.mat.valid.wt[, length(ss.mat.valid.wt[1, ])]) - 
                         min(ss.mat.valid.wt[, length(ss.mat.valid.wt[1, ])]))
    }
  } 
  
  # Valid Bhattacharyya distance results + Weight + Resampling probabilities.
  ss.mat.valid.wt.prob <- cbind(ss.mat.valid.wt,resamp.prob) 
  
  # Resample the indices. 
  resamp.ind <- sample(ss.mat.valid.wt.prob[, 1], size = length(paras[, 1]), 
                       replace = TRUE, prob = 
                       ss.mat.valid.wt.prob[, length(ss.mat.valid.wt.prob[1, ])]) 
  
  # Resampled parameter vectors, without perturbation.
  paras.nr.unperturbed <- paras[resamp.ind, ] 
  
  # An empty matrix used to store the perturbed parameter values.
  paras.nr.perturbed <- matrix(0, nrow = nrow(paras), ncol = ncol(paras)) 
  
  # Perturbation, see the manuscript for more details.
  for (i in 1:length(paras[1, ])) {
    for (j in 1:length(paras[, 1])){
      h = sqrt(1 - 0.05 ^ 2)
      paras.nr.perturbed[j, i] <- rnorm(1, h*paras.nr.unperturbed[j, i] + 
                                        (1 - h) * mean(paras.nr.unperturbed[, i]),
                                        0.05 * sd(paras.nr.unperturbed[, i]))
    }
  }  
  
  paras.nr.perturbed <- as.matrix(paras.nr.perturbed) 
  
  # Parameter values for the next round. 
  return(paras.nr.perturbed) 
  
}

ecm.rounds <- function(ref.paras, paras.init.ecm, n.rounds.ecm) {
# Purpose: Conduct the rounds of ECM-related trial.
# Inputs:
# ref.paras: Reference parameter values.
# paras.init.ecm: Initial parameter set of ECM-related trial.  
# n.rounds.ecm: Number of rounds being conducted.
# Output:
# bcd.res: Bhattacharyya distance results of all rounds.
# paras.res: Parameter sets of all rounds.

  n.row.paras <- length(paras.init.ecm[, 1])
  n.col.paras <- length(paras.init.ecm[1, ])
  
  bcd.mat.ecm <- matrix(0, nrow = n.rounds.ecm * n.row.paras, ncol = 2)
  paras.mat.ecm <- matrix(0, nrow = n.rounds.ecm * n.row.paras, ncol = n.col.paras)
  
  for (i in 1:n.rounds.ecm) {
    if (i == 1) {
      bcd.temp.ecm <- bcd(ref.paras, paras.init.ecm, ind = "ecm", 
                          print.bcd = TRUE)
      bcd.mat.ecm[1:n.row.paras, 1:2] <- bcd.temp.ecm
      
      paras.nr.ecm <- abc.bcd(bcd.temp.ecm, paras.init.ecm, rng.seed = 123)
      paras.mat.ecm[1:n.row.paras, 1:n.col.paras] <- paras.nr.ecm
    } else {
      bcd.temp.ecm <- bcd(ref.paras, paras.mat.ecm[((i - 2) * n.row.paras + 1):
                          ((i - 1) * n.row.paras), 1:n.col.paras], ind = "ecm", 
                          print.bcd = TRUE)
      bcd.mat.ecm[((i - 1) * n.row.paras + 1):(i * n.row.paras), 1:2] <- bcd.temp.ecm

      paras.nr.ecm <- abc.bcd(bcd.temp.ecm, paras.nr.ecm, rng.seed = 123)
      paras.mat.ecm[((i - 1) * n.row.paras + 1):(i * n.row.paras), 1:n.col.paras] <- paras.nr.ecm
    }
  }
  
  return(list(bcd.res = bcd.mat.ecm, paras.res = paras.mat.ecm))
}

ecm.mde.rounds <- function(ref.paras, paras.init.ecm.mde, n.rounds.ecm.mde) {
# Purpose: Conduct the rounds of ECM_MDE-related trial.
# Inputs:
# ref.paras: Reference parameter values.
# paras.init.ecm.mde: Initial parameter set of ECM_MDE-related trial.  
# n.rounds.ecm.mde: Number of rounds being conducted.
# Output:
# bcd.res: Bhattacharyya distance results of all rounds.
# paras.res: Parameter sets of all rounds.
  
  n.row.paras <- length(paras.init.ecm.mde[, 1])
  n.col.paras <- length(paras.init.ecm.mde[1, ])
  
  bcd.mat.ecm.mde <- matrix(0, nrow = n.rounds.ecm.mde * n.row.paras, ncol = 2)
  paras.mat.ecm.mde <- matrix(0, nrow = n.rounds.ecm.mde * n.row.paras, 
                              ncol = n.col.paras)
  
  for (i in 1:n.rounds.ecm.mde) {
    if (i == 1) {
      bcd.temp.ecm.mde <- bcd(ref.paras, paras.init.ecm.mde, ind = "ecm_mde", 
                              print.bcd = TRUE)
      bcd.mat.ecm.mde[1:n.row.paras, 1:2] <- bcd.temp.ecm.mde
      
      paras.nr.ecm.mde <- abc.bcd(bcd.temp.ecm.mde, paras.init.ecm.mde, 
                                  rng.seed = 123)
      paras.mat.ecm.mde[1:n.row.paras, 1:n.col.paras] <- paras.nr.ecm.mde
    } else {
      bcd.temp.ecm.mde <- bcd(ref.paras, paras.mat.ecm.mde[((i - 2) * 
                              n.row.paras + 1):((i - 1) * n.row.paras), 
                              1:n.col.paras], ind = "ecm_mde", 
                              print.bcd = TRUE)
      bcd.mat.ecm.mde[((i - 1) * n.row.paras + 1):(i * n.row.paras), 1:2] <- bcd.temp.ecm.mde
      
      paras.nr.ecm.mde <- abc.bcd(bcd.temp.ecm.mde, paras.nr.ecm.mde, 
                                  rng.seed = 123)
      paras.mat.ecm.mde[((i - 1) * n.row.paras + 1):(i * n.row.paras), 1:n.col.paras] <- paras.nr.ecm.mde
    }
  }
  
  return(list(bcd.res = bcd.mat.ecm.mde, paras.res = paras.mat.ecm.mde))
}

all3.rounds <- function(ref.paras, paras.init.all3, n.rounds.all3) {
# Purpose: Conduct the rounds of the final trial.
# Inputs:
# ref.paras: Reference parameter values.
# paras.init.all3: Initial parameter set of the final trial.  
# n.rounds.all3: Number of rounds being conducted.
# Output:
# bcd.res: Bhattacharyya distance results of all rounds.
# paras.res: Parameter sets of all rounds.
  
  n.row.paras <- length(paras.init.all3[, 1])
  n.col.paras <- length(paras.init.all3[1, ])
  
  bcd.mat.all3 <- matrix(0, nrow = n.rounds.all3 * n.row.paras, ncol = 2)
  paras.mat.all3 <- matrix(0, nrow = n.rounds.all3 * n.row.paras, 
                           ncol = n.col.paras)
  
  for (i in 1:n.rounds.all3) {
    if (i == 1) {
      bcd.temp.all3 <- bcd(ref.paras, paras.init.all3, ind = "all_3", 
                           print.bcd = TRUE)
      bcd.mat.all3[1:n.row.paras, 1:2] <- bcd.temp.all3
      
      paras.nr.all3 <- abc.bcd(bcd.temp.all3, paras.init.all3, 
                               rng.seed = 123)
      paras.mat.all3[1:n.row.paras, 1:n.col.paras] <- paras.nr.all3
    } else {
      bcd.temp.all3 <- bcd(ref.paras, paras.mat.all3[((i - 2) * n.row.paras + 
                           1):((i - 1) * n.row.paras), 1:n.col.paras], 
                           ind = "all_3", print.bcd = TRUE)
      bcd.mat.all3[((i - 1) * n.row.paras + 1):(i * n.row.paras), 1:2] <- bcd.temp.all3
      
      paras.nr.all3 <- abc.bcd(bcd.temp.all3, paras.nr.all3, rng.seed = 123)
      paras.mat.all3[((i - 1) * n.row.paras + 1):(i * n.row.paras), 1:n.col.paras] <- paras.nr.all3
    }
  }
  
  return(list(bcd.res = bcd.mat.all3, paras.res = paras.mat.all3))
}

