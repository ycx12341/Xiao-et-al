#### ABC-BCD scheme ########
### Authors: Yunchen Xiao & Len Thomas ###########

#### Set up ####
# Size of parameter sets
n.row.paras <- 10000

# Reference parameter values
ref.paras <- c(0.01, 0.05, 5, 10, 0.01, 0.1)

####### ECM initial paras #######################################
paras.ecm.r1 <- generate.init.paras(x = n.row.paras, rng.seed = 123)

write.results <- FALSE
if(write.results) {
  write.table(paras.ecm.r1, "Round 1 parameters 10000 ecm.txt")
}

###### ECM related rounds ######################################
n.rounds.ecm <- 3

ecm.results <- ecm.rounds(ref.paras, paras.init.ecm = paras.ecm.r1, 
                          n.rounds.ecm = n.rounds.ecm)

ecm.bcd.res <- ecm.results[["bcd.res"]]
write.results <- FALSE
if(write.results) {
  write.table(ecm.bcd.res, "ECM bcd all rounds.txt")
}

ecm.paras.res <- ecm.results[["paras.res"]]
write.results <- FALSE
if(write.results) {
  write.table(ecm.paras.res, "ECM paras all rounds.txt")
}
##### ECM + MDE initial paras #################################
ecm.paras.res.final <- ecm.paras.res[((n.rounds.ecm - 1) * n.row.paras + 1):
                       (n.rounds.ecm * n.row.paras), ]

# A different set.seed is used to avoid obtaining the same initial values for 
# the parameters which has not been covered by the scheme.
set.seed(121) 
dn <- runif(n.row.paras, 0.000069, 0.02)
gamma <- runif(n.row.paras, 0.005, 0.26)
eta <- ecm.paras.res.final[, 3]
dm <- runif(n.row.paras, 0.0001, 0.033)
alpha <- runif(n.row.paras, 0.07, 0.18)
rn <- runif(n.row.paras, 3.5, 9)

paras.ecm.mde.r1 <- cbind(dn, gamma, eta, dm, alpha, rn)
write.results <- FALSE
if(write.results) {
  write.table(paras.ecm.mde.r1, "Round 1 parameters 10000 ecm_mde.txt")
}

##### ECM + MDE related rounds ################################
n.rounds.ecm.mde <- 5
paras.init.ecm.mde <- paras.ecm.mde.r1

ecm.mde.results <- ecm.mde.rounds(ref.paras, 
                                  paras.init.ecm.mde = paras.ecm.mde.r1, 
                                  n.rounds.ecm.mde = n.rounds.ecm.mde)

ecm.mde.bcd.res <- ecm.mde.results[["bcd.res"]]
write.results <- FALSE
if(write.results) {
  write.table(ecm.mde.bcd.res, "ECM MDE bcd all rounds.txt")
}

ecm.mde.paras.res <- ecm.mde.results[["paras.res"]]
write.results <- FALSE
if(write.results) {
  write.table(ecm.mde.paras.res, "ECM MDE paras all rounds.txt")
}

##### All 3 initial paras #################################
ecm.mde.paras.res.final <- ecm.mde.paras.res[((n.rounds.ecm.mde - 1) * 
                           n.row.paras + 1):(n.rounds.ecm.mde * n.row.paras), ]

# A different set.seed is used to avoid obtaining the same initial parameters.
set.seed(120) 
dn <- runif(n.row.paras, 0.000069, 0.02)
gamma <- runif(n.row.paras, 0.005, 0.26)
eta <- ecm.mde.paras.res.final[, 3]
dm <- ecm.mde.paras.res.final[, 4]
alpha <- ecm.mde.paras.res.final[, 5]
rn <- runif(n.row.paras, 3.5, 9)

paras.all3.r1 <- cbind(dn, gamma, eta, dm, alpha, rn)
write.results <- FALSE
if(write.results) {
  write.table(paras.all3.r1, "Round 1 parameters 10000 all 3.txt")
}

##### All 3 related rounds ################################
n.rounds.all3 <- 9
all3.results <- all3.rounds(ref.paras, paras.init.all3 = paras.all3.r1, 
                            n.rounds.all3 = n.rounds.all3)

all3.bcd.res <- all3.results[["bcd.res"]]
write.results <- FALSE
if(write.results) {
  write.table(all3.bcd.res, "All 3 bcd all rounds.txt")
}

all3.paras.res <- all3.results[["paras.res"]]
write.results <- FALSE
if(write.results) {
  write.table(all3.paras.res, "All 3 paras all rounds.txt")
}
