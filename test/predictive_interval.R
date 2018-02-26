# Prediction Intervals
setwd('C:/Users/rbutl/OneDrive/Documents/NCBI Hachathon 2017/Clinotator/test/')
library(dplyr)
library(ggplot2)
library(rcompanion)
big_tbl <- read.csv('clinotator.clinvar_ge_2str.tsv', header=TRUE, sep='\t')

# US distribution
US_dist <- select(filter(big_tbl, grepl('Uncertain significance', CVCS) & CVNA >= 2),c(CTRS))
US_d2 <- sample_n(US_dist, 5000)
plotNormalHistogram(US_d2$CTRS, breaks = 50)
US_d2$CTRS = US_d2$CTRS - median(US_d2$CTRS)
plotNormalHistogram(US_d2$CTRS, breaks = 50)
US_d2$CTRS = sign(US_d2$CTRS) * (abs(US_d2$CTRS))^(1/10)
plotNormalHistogram(US_d2$CTRS, breaks = 50)
shapiro.test(US_d2$CTRS)
