# all_paras.R
# Author: Yunchen Xiao
# This R file reads in the parameter estimations at the end of different rounds
# of the ABC-BCD scheme.

paras_ori_ecm <- read.table("Round 1 parameters 10000 ecm.txt",sep="",header = TRUE)
paras_r2_ecm <- read.table("Round 2 parameters 10000 ecm.txt", sep = "", header = TRUE)
paras_r3_ecm <- read.table("Round 3 parameters 10000 ecm.txt", sep = "", header = TRUE)
paras_r4_ecm <- read.table("Round 4 parameters 10000 ecm.txt", sep = "", header = TRUE)

paras_ori_em <- read.table("Round 1 parameters 10000 ecm_mde.txt",sep="",header = TRUE)
paras_r2_em <- read.table("Round 2 parameters 10000 ecm_mde.txt",sep="",header = TRUE)
paras_r3_em <- read.table("Round 3 parameters 10000 ecm_mde.txt",sep="",header = TRUE)
paras_r4_em <- read.table("Round 4 parameters 10000 ecm_mde.txt",sep="",header = TRUE)
paras_r5_em <- read.table("Round 5 parameters 10000 ecm_mde.txt",sep="",header = TRUE)
paras_r6_em <- read.table("Round 6 parameters 10000 ecm_mde.txt",sep="",header = TRUE)

paras_ori_all3 <- read.table("Round 1 parameters 10000 all 3.txt",sep="",header = TRUE)
paras_r2_all3 <- read.table("Round 2 parameters 10000 all 3.txt",sep="",header = TRUE)
paras_r3_all3 <- read.table("Round 3 parameters 10000 all 3.txt",sep="",header = TRUE)
paras_r4_all3 <- read.table("Round 4 parameters 10000 all 3.txt",sep="",header = TRUE)
paras_r5_all3 <- read.table("Round 5 parameters 10000 all 3.txt",sep="",header = TRUE)
paras_r6_all3 <- read.table("Round 6 parameters 10000 all 3.txt",sep="",header = TRUE)
paras_r7_all3 <- read.table("Round 7 parameters 10000 all 3.txt",sep="",header = TRUE)
paras_r8_all3 <- read.table("Round 8 parameters 10000 all 3.txt",sep="",header = TRUE)
paras_r9_all3 <- read.table("Round 9 parameters 10000 all 3.txt",sep="",header = TRUE)
paras_r10_all3 <- read.table("Round 10 parameters 10000 all 3.txt",sep="",header = TRUE)
