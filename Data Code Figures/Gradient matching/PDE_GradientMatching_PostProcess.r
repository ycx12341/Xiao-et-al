# Post process the gradient matching scheme results, potentially from
# multiple runs

library(readr)

#in.dir is a vector with one row for each directory containing results to use
in.dir <- c("SimRes1")
#one colour for each directory
cols <- c("red")
n.runs <- length(in.dir)
cv <- c(0.01, 0.025, 0.05, 0.075, 0.10)
n.cvs <- length(cv)
true.values <- read_rds(paste0(in.dir[1], "/cv1_sim_res.rds"))$true.values
pars <- names(true.values)
n.pars <- length(pars)
mean.est <- perc.err <- mean.sd <- true.sd <- perc.err.sd <- 
  mean.boot <- boot.sd <- perc.err.boot.sd <-
  array(NA, dim = c(n.cvs, n.pars, n.runs), dimnames = list(cv, pars, 1:n.runs))
for(i in 1:n.cvs) {
  for(j in 1:n.runs) {
    res <- read_rds(paste0(in.dir[j], "/cv", i, "_sim_res.rds"))
    mean.est[i, , j] <- apply(res$par.ests, 2, mean)
    perc.err[i, , j] <- abs(mean.est[i, , j] - true.values) / true.values * 100
    mean.sd[i, , j] <- apply(res$sd.ests, 2, mean)
    true.sd[i, , j] <- apply(res$par.ests, 2, sd)
    perc.err.sd[i, , j] <- (mean.sd[i, , j] - true.sd[i, , j]) /
      true.sd[i, , j] * 100
    boot <- read_rds(paste0(in.dir[j], "/cv", i, "_boot_res.rds"))
    mean.boot[i, , j] <- apply(boot$boot.mean, 2, mean)
    boot.sd[i, , j] <- apply(boot$boot.sd, 2, mean)
    perc.err.boot.sd[i, , j] <- (boot.sd[i, , j] - true.sd[i, , j]) /
      true.sd[i, , j] * 100
  }
}

pdf("SimRes.pdf")
old.par <- par(no.readonly= TRUE)
par(mfrow = c(3,2), mar = c(4, 4, 4, 1) + 0.1)
# Plot of estimates vs true values
for(i in 1:n.pars){
  plot(cv, mean.est[, i, 1], ylab = paste0("Est (", pars[i], ")"), 
    ylim = range(c(mean.est[, i, ]), true.values[i]), type = "n", 
    main = pars[i])
  abline (h = true.values[i], lty = 2, col = "black", pch = 19)
  for(j in 1:n.runs)
    lines(cv, mean.est[, i, j], type = "b", col = cols[j], pch = 19)
}
# Plot of percentage error of estimates
for(i in 1:n.pars) {
  plot(cv, perc.err[, i, 1], ylab = paste0("abs % err (", pars[i], ")"), 
       ylim = c(0, max(perc.err)), type = "n", main = pars[i])
  for(j in 1:n.runs)
    lines(cv, perc.err[, i, j], type = "b", col = cols[j], pch = 19)
  
}
# Plot of estimated standard error vs true (estimated) values
for(i in 1:n.pars){
  plot(cv, mean.sd[, i, 1], ylab = paste0("SD (", pars[i], ")"), 
    ylim = range(c(mean.sd[, i, ], true.sd[, i, ])), type = "n", 
    main = pars[i])
  for(j in 1:n.runs) {
    lines(cv, mean.sd[, i, j], type = "b", col = cols[j])
    lines(cv, true.sd[, i, j], type = "b", col = cols[j], lty = 2)
  }
}
# Plot of percentage error in standard error estimates
for(i in 1:n.pars){
  plot(cv, perc.err.sd[, i, 1], ylab = paste0("% err SD (", pars[i], ")"), 
       ylim = range(perc.err.sd), type = "n", main = pars[i])
  for(j in 1:n.runs)
    lines(cv, perc.err.sd[, i, j], type = "b", col = cols[j])
}
# Plot of estimated standard error from bootstrap vs true (estimated) values
for(i in 1:n.pars){
  plot(cv, boot.sd[, i, 1], ylab = paste0("bootSD (", pars[i], ")"), 
       ylim = range(c(boot.sd[, i, ], true.sd[, i, ])), type = "n", 
       main = pars[i])
  for(j in 1:n.runs) {
    lines(cv, boot.sd[, i, j], type = "b", col = cols[j])
    lines(cv, true.sd[, i, j], type = "b", col = cols[j], lty = 2)
  }
}
# Plot of percentage error in bootstrap standard error estimates
for(i in 1:n.pars){
  plot(cv, perc.err.boot.sd[, i, 1], ylab = paste0("% err bootSD (", pars[i], ")"), 
       ylim = range(perc.err.boot.sd), type = "n", main = pars[i])
  abline(h = 0, lty = 2)
  for(j in 1:n.runs)
    lines(cv, perc.err.boot.sd[, i, j], type = "b", col = cols[j])
}
# Plot of bootstrap mean vs true value
for(i in 1:n.pars){
  plot(cv, mean.boot[, i, 1], ylab = paste0("bootEst (", pars[i], ")"), 
       ylim = range(c(mean.boot[, i, ]), true.values[i]), type = "n", 
       main = pars[i])
  abline (h = true.values[i], lty = 2, col = "black", pch = 19)
  for(j in 1:n.runs)
    lines(cv, mean.boot[, i, j], type = "b", col = cols[j], pch = 19)
}

par <- old.par
dev.off()

