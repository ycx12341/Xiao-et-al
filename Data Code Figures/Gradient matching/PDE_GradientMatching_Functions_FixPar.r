#### Gradient matching scheme - Functions ########
### Author: Yunchen Xiao & Len Thomas ###########

generate.reference.data <- function(n.x11, max.x11, dt, max.t,
  dn, gamma, ita, dm, alpha, r, beta, n0, f0, m0) {
#Purpose: Generates a realization from the cancer invasion model
# Saves output at n.x11 space points evenly distributed between 0 and max.x11
# and at integer seconds from 0 to max.t.  Model run at time discretization dt
#Inputs:
# n.x11 - number of space points to evaluate model and save at
# max.x11 - max space value (minimum assumed to be 0)
# dt - time interval for discretization
# max.t - max time values (minimum assumed to be 0)
# dn, gamma, ita, dm, alpha, r, beta - model parameters
# n0, f0, m0 - vectors of initial conditions (of length n.x11)
#Outputs:
# n - tumor cell density
#     matrix containing values of n at (max.t + 1) x n.11 time and space points
# f - extracelluar matrix (ECM)
#     matrix containing values of f at (max.t + 1) x n.11 time and space points
# m - matrix degrading enzyme (ECM)
#     matrix containing values of m at (max.t + 1) x n.11 time and space points
#Implementation note:
# At present restricted to return data from integer time points.  Would be
# worthwhile to extend to return at arbitrary time points (within the level
# of precision set by dt)
  
  #Calculate space points and space discretization
  x11 <- seq(0, max.x11, length = n.x11)
  h <- 1/(n.x11 - 1)
  #Calculate time interval to save results at
  t.to.save <- round(1/dt)

  # Store the initial values into the result matrix
  res.n <- res.f <- res.m <- matrix(0, max.t + 1, n.x11)
  res.n [1, ] <- n0
  res.f [1, ] <- f0
  res.m [1, ] <- m0

  # Run the equations
  p <- 1
  t <- 2
  while(p * dt <= max.t) {
    f[2:(n.x11 - 1)] <- -ita * dt * m[2:(n.x11 - 1)] * f[2:(n.x11 - 1)] + 
      f[2:(n.x11 - 1)]

    m[2:(n.x11 - 1)] <- dm * (m[1:(n.x11 - 2)] + m[3:n.x11] - 2 * 
                                m[2:(n.x11 - 1)]) * dt / (h ^ 2) +
      alpha * n[2:(n.x11 - 1)] * dt -  
      beta * m[2:(n.x11 - 1)] * dt + m[2:(n.x11 - 1)]
    
    n[2:(n.x11 - 1)] <- dn * (n[1:(n.x11 - 2)] + n[3:n.x11] - 2 * 
                                n[2:(n.x11 - 1)]) * dt / (h ^ 2) -  
      gamma * (n[3:n.x11] - n[2:(n.x11 - 1)]) * 
        (f[3:n.x11] - f[2:(n.x11 - 1)]) * dt / (h ^ 2) - 
      gamma * n[2:(n.x11 - 1)] * (f[1:(n.x11 - 2)] + f[3:n.x11] - 
                                    2 * f[2:(n.x11 - 1)]) * dt / (h ^ 2) + 
      r * (1 - f[2:(n.x11 - 1)] - n[2:(n.x11 - 1)]) * n[2:(n.x11 - 1)] * dt + 
        n[2:(n.x11 - 1)]
  
    #No flux boundary condition
    n[1] <- n[2]
    n[n.x11] <- n[n.x11 - 1]
    
    f[1] <- f[2]
    f[n.x11] <- f[n.x11 - 1]
    
    m[1] <- m[2]
    m[n.x11] <- m[n.x11 - 1]
  
    # Save the results at integer time steps
    if(p %% t.to.save == 0) {
      res.n[t, ] <- n; res.f[t, ] <- f; res.m[t, ] <- m
      t <- t + 1
    }
  
    p <- p + 1
  }
  return(list (n = res.n, f = res.f, m = res.m))

}

perturb.reference.data <- function(ref.data, cv = NULL, sd = NULL, distribution = "gamma") {
#Purpose: Adds random error to the input dataset
# Note - for input values less than or equal to zero, no error is added
#Inputs:
# ref.data - list with 3 elements n, f and m, each time x space matrices
# cv - desired coefficient of variation - must be specified for distribution = "gamma"
#   and either this or sd should be specified for distribution = "normal"
# sd - desired standard deviation - an alternative to specifying cv for normal distribution
#   data
# distribution - if "gamma" adds gamma error; if "normal" adds normal error
#Outputs:
# named list of same structure as ref.data, containing perturbed values
  
  #Define function that adds gamma or normal error to a vector of values
  perturb.vec <- function(vec, cv, distribution) {
    ind <- (vec > 0)
    n <- length(vec[ind])
    if(distribution == "gamma") {
      if(is.null(cv)) stop("cv must be specified for gamma distribution\n")
      shape <- (1 / cv)^2
      vec[ind] <- rgamma(n, shape, rate = shape / vec[ind])
    } else {
      if(distribution == "normal") {
        if(is.null(cv)&is.null(sd)) 
          stop("cv or sd must be specified for normal distribution\n")
        if(!is.null(cv)) {
          vec[ind] <- rnorm(n, vec[ind], vec[ind] * cv)
        } else {
          vec[ind] <- rnorm(n, vec[ind], sd)
        }
      } else {
        stop("Specified distribution is not normal or gamma\n")
      }
    }
    return(vec)
  }

  #Apply the perturbation
  ref.data$n <- apply(ref.data$n, 2, perturb.vec, cv, distribution)
  ref.data$f <- apply(ref.data$f, 2, perturb.vec, cv, distribution)
  ref.data$m <- apply(ref.data$m, 2, perturb.vec, cv, distribution)
  return(ref.data)
}

approximate.gradients <- function(data, x11, max.t, distribution = "gamma") {
  #Purpose: Returns several first and second derivatives of the 3 input time series
  #Inputs:
  # data - list containing elements
  #   n - matrix containing values of n at time and space points
  #   f - matrix containing values of f at time and space points
  #   m - matrix containing values of m at time and space points
  # x11 - vector of space points
  # max.t - max time point (time points assumed to start with 0 and be integer)
  # distribution - if "gamma" then assumes gamma family for the gam, with log link;
  #  if "normal" assumes normal family for gam, with identity link
  #Outputs:
  # grad.lhs_n: Temporal gradients of n. (Tumour cells)
  # grad.rhs_dn: Spatial gradients associated with dn. (2nd order derivative of n in space.)
  # grad.rhs_gamma: Spatial gradients associated with gamma. (Chemotaxis terms.)
  # grad.rhs_r: Spatial gradients associated with rn. (Density of tumour cells, n.)
  # grad.lhs_f: Temporal gradients of f. (ECM)
  # grad.rhs_ita: Spatial gradients associated with ita. (Product of ECM and MDE densities.)
  # grad.lhs_m: Temporal gradients of m. (MDE)
  # grad.rhs_dm: Spatial gradients associated with dm. (2nd order derivative of m in space.)
  # grad.rhs_alpha: Spatial gradients associated with alpha. (Density of tumour cells, n.)
  
  #mgcv package required for fitting gams
  require(mgcv) 
  
  #Space discretization for calculating derivatives is set by data spacing 
  # (but need not be)
  h <- x11[2] - x11[1]
  #Time discretization for calculating derivatives is set manually
  dt <- 0.001
  
  #Transform each time series into a data frame
  n.x11 <- length(x11)
  tp_trun <- seq(1, max.t - 1, by = 1)
  dat_n <- data.frame(n = as.vector(data$n), t = rep(tp_trun, times = n.x11 - 2), 
                      x11 = rep(x11[2:(n.x11 - 1)], each = max.t - 1))
  dat_f <- data.frame(f = as.vector(data$f), t = rep(tp_trun, times = n.x11 - 2), 
                      x11 = rep(x11[2:(n.x11 - 1)], each = max.t - 1))
  dat_m <- data.frame(m = as.vector(data$m), t = rep(tp_trun, times = n.x11 - 2), 
                      x11 = rep(x11[2:(n.x11 - 1)], each = max.t - 1))
  
  #Fit a gam to each time series
  if(distribution == "gamma") {
    fam <- Gamma(link = "log")
  } else {
    if(distribution == "normal") {
      fam <- gaussian(link = "identity")
    } else {
      stop("Distribution not one of gamma or normal\n")
    }
  }
  suppressWarnings(spl <- gam(n ~ s(t, x11, bs = "ad"), family = fam, data = dat_n))
  suppressWarnings(spl2 <- gam(f ~ s(t, x11, bs = "ad"), family = fam, data = dat_f))
  suppressWarnings(spl3 <- gam(m ~ s(t, x11, bs = "ad"), family = fam, data = dat_m))
  
  #Prepare to calculate the required gradients
  tp <- seq(0, max.t, by = 1)
  dim_row <- length(data$n[, 1])
  dim_col <- length(data$f[1, ])
  
  # Temporal gradients of ECM.
  grad_lhs_f_data <- matrix(0, dim_row, dim_col)
  for (i in 1:dim_row) {
    for (j in 1:dim_col) {
      predict_f_t1 <- predict(spl2, newdata = data.frame(t = tp[i + 1] + dt, x11 = x11[j + 1]), 
                              type = "response")
      predict_f_t2 <- predict(spl2, newdata = data.frame(t = tp[i + 1] - dt, x11 = x11[j + 1]), 
                              type = "response")
      grad_lhs_f_data[i,j] <- (predict_f_t1 - predict_f_t2) / (2 * dt)
    }
  }
  
  # Spatial gradients of ECM.
  grad_rhs_f_ita <- matrix(0, dim_row, dim_col)
  for (i in 1:dim_row) {
    for (j in 1:dim_col) {
      predict_f <- predict(spl2, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]), 
                           type = "response")
      predict_m <- predict(spl3, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]), 
                           type = "response")
      grad_rhs_f_ita[i, j] <- predict_f * predict_m
    }
  }
  
  # Temporal gradients of MDE.  
  grad_lhs_m_data <- matrix(0, dim_row, dim_col)
  for (i in 1:dim_row) {
    for (j in 1:dim_col) {
      predict_m_t1 <- predict(spl3, newdata = data.frame(t = tp[i + 1] + dt, x11 = x11[j + 1]), 
                              type = "response")
      predict_m_t2 <- predict(spl3, newdata = data.frame(t = tp[i + 1] - dt, x11 = x11[j + 1]), 
                              type = "response")
      grad_lhs_m_data[i,j] <- (predict_m_t1 - predict_m_t2) / (2 * dt)
    }
  }
  
  # Spatial gradients of MDE.
  grad_rhs_m_dm <- matrix(0, dim_row, dim_col)
  for (i in 1:dim_row) {
    for (j in 1:dim_col) {
      predict_x1 <- predict(spl3, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] - h), 
                            type = "response")
      predict_x2 <- predict(spl3, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]), 
                            type = "response")
      predict_x3 <- predict(spl3, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] + h), 
                            type = "response")
      grad_rhs_m_dm[i, j] <- (predict_x3 + predict_x1 - 2 * predict_x2) / (h ^ 2)
    }
  }
  grad_rhs_m_alpha <- matrix(0, dim_row, dim_col)
  for (i in 1:dim_row) {
    for (j in 1:dim_col) {
      predict_n <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]), 
                           type = "response")
      grad_rhs_m_alpha[i, j] <- predict_n
    }
  }
  
  # Temporal gradients of tumour cells.
  grad_lhs_n_data <- matrix(0, dim_row, dim_col)
  for (i in 1:dim_row) {
    for (j in 1:dim_col) {
      predict_n_t1 <- predict(spl, newdata = data.frame(t = tp[i + 1] + dt, x11 = x11[j + 1]), 
                              type = "response")
      predict_n_t2 <- predict(spl, newdata = data.frame(t = tp[i + 1] - dt, x11 = x11[j + 1]), 
                              type = "response")
      grad_lhs_n_data[i,j] <- (predict_n_t1 - predict_n_t2) / (2 * dt)
    }
  }
  
  # Spatial gradients of tumour cells.
  grad_rhs_n_dn <- matrix(0, dim_row, dim_col)
  for (i in 1:dim_row) {
    for (j in 1:dim_col) {
      predict_x1 <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] - h), 
                            type = "response")
      predict_x2 <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]), 
                            type = "response")
      predict_x3 <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] + h), 
                            type = "response")
      grad_rhs_n_dn[i, j] <- (predict_x3 + predict_x1 - 2 * predict_x2) / (h ^ 2)
    }
  }
  grad_rhs_n_gamma <- matrix(0, dim_row, dim_col)
  for (i in 1:dim_row) {
    for (j in 1:dim_col) {
      predict_x1_n <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] - h), 
                              type = "response")
      predict_x2_n <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]), 
                              type = "response")
      predict_x3_n <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] + h), 
                              type = "response")
      predict_n.dash <- (predict_x3_n - predict_x1_n) / (2 * h)
      predict_x1_f <- predict(spl2, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] - h), 
                              type = "response")
      predict_x2_f <- predict(spl2, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]), 
                              type = "response")
      predict_x3_f <- predict(spl2, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] + h), 
                              type = "response")
      predict_f.dash <- (predict_x3_f - predict_x1_f) / (2 * h)
      predict_f.ddash <- (predict_x3_f + predict_x1_f - 2 * predict_x2_f) / (h ^ 2)
      grad_rhs_n_gamma[i, j] <- predict_n.dash * predict_f.dash + predict_x2_n * predict_f.ddash
    }
  }
  grad_rhs_n_r <- matrix(0, dim_row, dim_col)
  for (i in 1:dim_row) {
    for (j in 1:dim_col) {
      predict_n <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]),
                           type = "response")
      predict_f <- predict(spl2, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]),
                           type = "response")
      grad_rhs_n_r[i, j] <- predict_n * (1 - predict_n - predict_f)
    }
  }
  
  #Return required gradients as a named list  
  return(list(grad.lhs_n = grad_lhs_n_data, grad.rhs_dn = grad_rhs_n_dn, 
              grad.rhs_gamma = grad_rhs_n_gamma, grad.rhs_r = grad_rhs_n_r, 
              grad.lhs_f = grad_lhs_f_data, grad.rhs_ita = grad_rhs_f_ita,
              grad.lhs_m = grad_lhs_m_data, grad.rhs_dm = grad_rhs_m_dm, 
              grad.rhs_alpha = grad_rhs_m_alpha))
}

calculate.sse <- function(par, grads, fixed.par = rep(NA, 6), debug = FALSE) {
  #Purpose: Return the objective function -- sum of squared errors (differences) 
  #between lhs and rhs of model equations, given approximated gradients
  #Inputs:
  # par - model parameter values being estimated, in a vector.  Parameters
  #  that are part of the model, but not being estimated, are omitted
  # grads - list of model gradients (see approximate.gradients outputs for elements)
  # fixed.par - vector of length 6 corresponding to all parameters in the model.
  #  Parameters that are estimated are NA, while those that are fixed have their
  #  values here.  The number of NA entries should therefore correspond with the 
  #  length of par.
  # debug - if TRUE, prints out the objective function
  #Outputs:
  # value of objective function
  
  #Combine parameters to estimate and those that are fixed into one vector
  fixed.par[is.na(fixed.par)] <- par
  
  #Calculate SSE for each equation     
  sse1 <- sum((grads$grad.lhs_n - (fixed.par[1] * grads$grad.rhs_dn - 
                                     fixed.par[2] * grads$grad.rhs_gamma + par[3] * grads$grad.rhs_r))^2)
  sse2 <- sum((grads$grad.lhs_f + fixed.par[4] * grads$grad.rhs_ita)^2)
  sse3 <- sum((grads$grad.lhs_m - (fixed.par[5] * grads$grad.rhs_dm + 
                                     fixed.par[6] * grads$grad.rhs_alpha))^2)
  
  # Sum the sse from each equation
  res <- sum(sse1, sse2, sse3)
  if (debug) cat(res, par, "\n")
  return(res)
}

estimate.sd <- function (optim.res, df) {
  #Purpose: Returns standard error of parameters from a least-squares fit
  # using optim
  #Inputs:
  # optim.res - output from calling optim (with option hessian = TRUE)
  # df - degrees of freedom for fit (number of data points - number of params)
  #Implementation note: if hessian can not be inverted, returns vector of NAs for SDs
  
  resid.var <- optim.res$value / df
  sd.ests <- try(sqrt(resid.var * 2 * diag(solve(optim.res$hessian))), silent = TRUE)
  if(class(sd.ests) == "try-error") sd.ests <- rep(NA, length(optim.res$par))
  return(sd.ests)
}
