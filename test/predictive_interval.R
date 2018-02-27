# Prediction Intervals (Non-parametric)
setwd('C:/Users/rbutl/OneDrive/Documents/NCBI Hachathon 2017/Clinotator/test/')
library(dplyr)
library(ggplot2)
big_tbl <- read.csv('clinotator.clinvar_ge_2str.tsv', header=TRUE, sep='\t')

# US distribution
US_tbl <- filter(big_tbl, grepl('Uncertain significance', CVCS) & CVNA >= 2)
US_dist <- sort(US_tbl$CTRS)
US_cth <- floor((0.01/2) * (length(US_dist) + 1))
US_PI <- c(US_dist[US_cth], US_dist[length(US_dist) - US_cth + 1])

# LP distribution
LP_tbl <- filter(big_tbl, grepl('^Likely pathogenic', CVCS) & CVNA >= 2)
LP_dist <- sort(LP_tbl$CTRS)
LP_cth <- floor((0.01/2) * (length(LP_dist) + 1))
LP_PI <- c(LP_dist[LP_cth], LP_dist[length(LP_dist) - LP_cth + 1])

# LB distribution
LB_tbl <- filter(big_tbl, grepl('^Likely benign', CVCS) & CVNA >= 2)
LB_dist <- sort(LB_tbl$CTRS)
LB_cth <- floor((0.01/2) * (length(LB_dist) + 1))
LB_PI <- c(LB_dist[LB_cth], LB_dist[length(LB_dist) - LB_cth + 1])

# P distribution
P_tbl <- filter(big_tbl, grepl('Pathogenic($|,)', CVCS) & CVNA >= 2)
P_dist <- sort(P_tbl$CTRS)
P_cth <- floor((0.01/2) * (length(P_dist) + 1))
P_PI <- c(P_dist[P_cth], P_dist[length(P_dist) - P_cth + 1])

# B distribution
B_tbl <- filter(big_tbl, grepl('Benign($|,)', CVCS) & CVNA >= 2)
B_dist <- sort(B_tbl$CTRS)
B_cth <- floor((0.01/2) * (length(B_dist) + 1))
B_PI <- c(B_dist[B_cth], B_dist[length(B_dist) - B_cth + 1])

cat('B:', median(B_dist), B_PI, '\n')
cat('LB:', median(LB_dist), LB_PI, '\n')
cat('US:', median(US_dist), US_PI, '\n')
cat('LP:', median(LP_dist), LP_PI, '\n')
cat('P:', median(P_dist), P_PI, '\n')

