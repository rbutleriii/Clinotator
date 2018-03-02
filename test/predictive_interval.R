# Prediction Intervals (Non-parametric)
setwd('C:/Users/rbutl/OneDrive/Documents/NCBI Hachathon 2017/Clinotator/test/')
library(dplyr)
library(ggplot2)
big_tbl <- read.csv('clinotator.clinvar_ge_2str.tsv', header=TRUE, sep='\t', na.strings='.')

# US distribution
US_tbl <- filter(big_tbl, grepl('Uncertain significance', CVCS) & CVNA >= 2 & CVSZ == 2)
US_dist <- sort(US_tbl$CTRS)
US_cth <- floor((0.01/2) * (length(US_dist) + 1))
US_PI <- c(US_dist[US_cth], US_dist[length(US_dist) - US_cth + 1])
US_quant <- quantile(US_dist, probs = c(0.01, 0.99))

# LP distribution
LP_tbl <- filter(big_tbl, grepl('^Likely pathogenic', CVCS) & CVNA >= 2 & CVSZ == 2)
LP_dist <- sort(LP_tbl$CTRS)
LP_cth <- floor((0.03/2) * (length(LP_dist) + 1))
LP_PI <- c(LP_dist[LP_cth], LP_dist[length(LP_dist) - LP_cth + 1])
LP_quant <- quantile(LP_dist, probs = c(0.02, 0.98))

# LB distribution
LB_tbl <- filter(big_tbl, grepl('^Likely benign', CVCS) & CVNA >= 2 & CVSZ == 2)
LB_dist <- sort(LB_tbl$CTRS)
LB_cth <- floor((0.005/2) * (length(LB_dist) + 1))
LB_PI <- c(LB_dist[LB_cth], LB_dist[length(LB_dist) - LB_cth + 1])
LB_quant <- quantile(LB_dist, probs = c(0.003, 0.997))

# P distribution
P_tbl <- filter(big_tbl, grepl('Pathogenic($|,)', CVCS) & CVNA >= 2 & CVSZ == 2)
P_dist <- sort(P_tbl$CTRS)
P_cth <- floor((0.04/2) * (length(P_dist) + 1))
P_PI <- c(P_dist[P_cth], P_dist[length(P_dist) - P_cth + 1])
P_quant <- quantile(P_dist, probs = c(0.02, 0.98))

# B distribution
B_tbl <- filter(big_tbl, grepl('Benign($|,)', CVCS) & CVNA >= 2 & CVSZ == 2)
B_dist <- sort(B_tbl$CTRS)
B_cth <- floor((0.05/2) * (length(B_dist) + 1))
B_PI <- c(B_dist[B_cth], B_dist[length(B_dist) - B_cth + 1])
B_quant <- quantile(B_dist, probs = c(0.03, 0.97))

cat('B:', median(B_dist), B_PI, 'Quant:', B_quant, '\n')
cat('LB:', median(LB_dist), LB_PI, 'Quant:', LB_quant, '\n')
cat('US:', median(US_dist), US_PI, 'Quant:', US_quant, '\n')
cat('LP:', median(LP_dist), LP_PI, 'Quant:', LP_quant, '\n')
cat('P:', median(P_dist), P_PI, 'Quant:', P_quant, '\n')

##Calibration Values
# -6, -3, -0.3, 3, 6
# B: -12 -39 -8.4 Quant: -37.2 -8.4 
# LB: -6 -24 -4.2 Quant: -23.7 -4.2 
# US: -0.6 -2.37 -0.42 Quant: -2.1 -0.45 
# LP: 6 4.2 15 Quant: 4.2 12.018 
# P: 12 8.4 37.2 Quant: 8.4 37.2 

## Initial Testing 
# -5, -3, -0.3, 1.6, 2.9
# B: -10 -43.5 -6 Quant: -38 -6.5 
# LB: -6 -21 -4.5 Quant: -17.7 -4.5 
# US: -0.6 -2.37 -0.42 Quant: -2.1 -0.45 
# LP: 3.2 1.6 10.08 Quant: 2.08 8 
# P: 5.8 2.9 25.23 Quant: 3.48 22.04 

# -5.5, -2.9, -0.3, 2, 4
# B: -11 -47.85 -6.6 Quant: -41.8 -7.15 
# LB: -5.8 -20.3 -4.35 Quant: -17.11 -4.35 
# US: -0.6 -2.37 -0.42 Quant: -2.1 -0.45 
# LP: 4 2 12.6 Quant: 2.6 10 
# P: 8 4 34.8 Quant: 4.8 30.4 

# -6.5, -3.9, -0.3, 3, 5
# B: -13 -56.55 -7.8 Quant: -49.4 -8.45 
# LB: -7.8 -27.3 -5.85 Quant: -23.01 -5.85 
# US: -0.6 -2.37 -0.42 Quant: -2.1 -0.45 
# LP: 6 3 18.9 Quant: 3.9 15 
# P: 10 5 43.5 Quant: 6 38 

# -8, -5, -0.3, 6, 8
# B: -16 -69.6 -9.6 Quant: -60.8 -10.4 
# LB: -10 -35 -7.5 Quant: -29.5 -7.5 
# US: -0.6 -2.37 -0.42 Quant: -2.1 -0.45 
# LP: 12 6 37.8 Quant: 7.8 30 
# P: 16 8 69.6 Quant: 9.6 60.8 
