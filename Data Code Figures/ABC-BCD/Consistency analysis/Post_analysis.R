# Post_analysis.R
# Author: Yunchen Xiao
# This R file calculates the Monte-Carlo errors and error percentage of the
# parameter estimations obtained from the three different runs of applying ABC-BCD scheme
# on the main reference dataset.   

# Read in the final parameter estimations from each run
paras_post_ori <- read.table("Round 10 parameters 10000 all 3_ori.txt",sep="",header = TRUE)
paras_post_ex1 <- read.table("Round 10 parameters 10000 all 3_ex1.txt",sep="",header = TRUE)
paras_post_ex2 <- read.table("Round 10 parameters 10000 all 3_ex2.txt",sep="",header = TRUE)

dn_final <- c(mean(paras_post_ori[,1]), mean(paras_post_ex1[,1]), mean(paras_post_ex2[,1]))
ga_final <- c(mean(paras_post_ori[,2]), mean(paras_post_ex1[,2]), mean(paras_post_ex2[,2]))
eta_final <- c(mean(paras_post_ori[,3]), mean(paras_post_ex1[,3]), mean(paras_post_ex2[,3]))
dm_final <- c(mean(paras_post_ori[,4]), mean(paras_post_ex1[,4]), mean(paras_post_ex2[,4]))
al_final <- c(mean(paras_post_ori[,5]), mean(paras_post_ex1[,5]), mean(paras_post_ex2[,5]))
rn_final <- c(mean(paras_post_ori[,6]), mean(paras_post_ex1[,6]), mean(paras_post_ex2[,6]))

# Monte-Carlo errors of each parameter 
sd_dn <- sd(dn_final)/sqrt(3) # 0.0001439805
sd_ga <- sd(ga_final)/sqrt(3) # 0.0006608775
sd_eta <- sd(eta_final)/sqrt(3) # 0.05611386
sd_dm <- sd(dm_final)/sqrt(3) # 2.318548e-05
sd_al <- sd(al_final)/sqrt(3) # 0.0002189775
sd_rn <- sd(rn_final)/sqrt(3) # 0.06278749

# Corresponding error percentages.
sd(dn_final)/0.01*100 # 2.49%
sd(ga_final)/0.05*100 # 2.29%
sd(rn_final)/5*100 # 2.18%
sd(eta_final)/10*100 # 0.972%
sd(dm_final)/0.01*100 # 0.402%
sd(al_final)/0.1*100 # 0.379%
