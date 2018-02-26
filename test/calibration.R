# protocol for BCa resampling to calibrate Clinotator
setwd('C:/Users/rbutl/OneDrive/Documents/NCBI Hachathon 2017/Clinotator/test/')
library(dplyr)
library(ggplot2)
library(simpleboot)
library(nortest)
big_tbl <- read.csv('clinotator.clinvar_ge_2str.tsv', header=TRUE, sep='\t')

# P Set
P_dist <- select(filter(big_tbl, grepl('Pathogenic($|,)', CVCS) & CVNA >= 2),c(CTRS))
P_dist_n <- length(P_dist$CTRS)
ggplot(P_dist, aes(CTRS)) + geom_density()
P_boot <- one.boot(P_dist$CTRS, mean, R=(P_dist_n + 1), tr=.10)
P_boot_dist <- data.frame(P_boot$t)
ggplot(P_boot_dist, aes(x=P_boot.t)) + geom_density() + geom_vline(aes(xintercept=P_boot$t0), col='red')
ad.test(P_boot_dist$P_boot.t)
boot.ci(P_boot, type=c("bca"))

# LP Set
LP_dist <- select(filter(big_tbl, grepl('Pathogenic($|,)', CVCS) & CVNA >= 2),c(CTRS))
LP_dist_n <- length(LP_dist$CTRS)
ggplot(LP_dist, aes(CTRS)) + geom_density()
LP_boot <- one.boot(LP_dist$CTRS, mean, R=(LP_dist_n + 1), tr=.10)
LP_boot_dist <- data.frame(LP_boot$t)
ggplot(LP_boot_dist, aes(x=LP_boot.t)) + geom_density() + geom_vline(aes(xintercept=LP_boot$t0), col='red')
ad.test(LP_boot_dist$LP_boot.t)
boot.ci(LP_boot, type=c("bca"))

# US Set
US_dist <- select(filter(big_tbl, grepl('Uncertain significance', CVCS) & CVNA >= 2), c(CTRS))
US_dist_n <- length(US_dist$CTRS)
ggplot(US_dist, aes(CTRS)) + geom_density()
US_boot <- one.boot(US_dist$CTRS, mean, R=(US_dist_n + 1), tr=.10)
US_boot_dist <- data.frame(US_boot$t)
ggplot(US_boot_dist, aes(x=US_boot.t)) + geom_density() + geom_vline(aes(xintercept=US_boot$t0), col='red')
ad.test(US_boot_dist$US_boot.t)
boot.ci(US_boot, type=c("bca"), conf = 0.99)

# LB Set
LB_dist <- select(filter(big_tbl, grepl('Pathogenic($|,)', CVCS) & CVNA >= 2),c(CTRS))
LB_dist_n <- length(LB_dist$CTRS)
ggplot(LB_dist, aes(CTRS)) + geom_density()
LB_boot <- one.boot(LB_dist$CTRS, mean, R=(LB_dist_n + 1), tr=.10)
LB_boot_dist <- data.frame(LB_boot$t)
ggplot(LB_boot_dist, aes(x=LB_boot.t)) + geom_density() + geom_vline(aes(xintercept=LB_boot$t0), col='red')
ad.test(LB_boot_dist$LB_boot.t)
boot.ci(LB_boot, type=c("bca"))

# B Set
B_dist <- select(filter(big_tbl, grepl('Pathogenic($|,)', CVCS) & CVNA >= 2),c(CTRS))
B_dist_n <- length(B_dist$CTRS)
ggplot(B_dist, aes(CTRS)) + geom_density()
B_boot <- one.boot(B_dist$CTRS, mean, R=(B_dist_n + 1), tr=.10)
B_boot_dist <- data.frame(B_boot$t)
ggplot(B_boot_dist, aes(x=B_boot.t)) + geom_density() + geom_vline(aes(xintercept=B_boot$t0), col='red')
ad.test(B_boot_dist$B_boot.t)
boot.ci(B_boot, type=c("bca"))