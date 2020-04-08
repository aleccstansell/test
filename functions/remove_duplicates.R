#####################################################################################################################################
#Remove duplicates function - removes duplicate ensembl gene ids produced from Deseq2
#####################################################################################################################################


#Load dplyr
library(dplyr)
#####################################################################
#Function
#####################################################################

removedup <- function(dataframe){
  
  nodup <- distinct(dataframe)
  
  return(nodup)
}
######################################################################
# setwd("~/doorknobdave/alecS/results/maydeseq/pairwise_deseq_with_biomart")
# 
# data <- read.csv("NvsCM.csv")
# 
# rem <- unique(data)
# 
# write.csv(rem, "NvsCM_dup_rem_sig.csv")
# 
# 
# 
# 
# ?unique
# 
# data
# 
# unique <- unique(data)
# 
# write.csv(unique, "nodupNvsCMresbiomart1.csv")
# return(nodup)
# 
# read.csv()
# 
# 
# 
# 
# 
