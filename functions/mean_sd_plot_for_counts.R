#
'Code to plot mean sd plot for different normalised coutnts'
#

#load packages
BiocManager::install("vsn")
library(vsn)

#Load Normalised counts
setwd("~/doorknobdave/alecS/alec_code/coexpression/final_coexpression/data")
load("final_heat_normalised.RData")

#Plot meansSD plots

options(repr.plot.height=3, repr.plot.width=5)
log.cts.one <- as.matrix(log2(normalised_clean$input_counts))


#log.cts.one <- log.cts.one[rowSums(log.cts.one > 1), ]
siszefact <- meanSdPlot(as.matrix(log2(normalised_clean$size_factor_counts), ranks = FALSE, bins = 120))
dev.off()
varstab <- meanSdPlot(as.matrix(log2(normalised_clean$variance_stabilised_counts)), ranks = FALSE, bins = 120)
varstab$gg + ylim(0, 5)
raw <- meanSdPlot((abs(as.matrix(log2(normalised_clean$input_counts)))), ranks = FALSE, bins = 120)
raw$gg + ylim(0,5) #+ xlim(0, 17)

rlog <- meanSdPlot(as.matrix(log2(normalised_clean$size_factor_counts)), ranks = FALSE, bins = 100)

rlog$gg 


library("ggplot2")
siszefact$gg #+ xlim(1, 4)
varstab$gg 

raw$gg + xlim(0, 16)
rlog$gg + ylim(0, 5)

gg



?meanSdPlot
dev.off()
