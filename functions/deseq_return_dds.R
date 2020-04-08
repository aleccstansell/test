####################################################################################
#Deseq With Graphs
#Performs Deseq for a single comparion that is chosen by the user
####################################################################################


#Load Libraries
library(praise)
library(tximport)
library(DESeq2)
library(tidyverse)
library(rhdf5)

deseq_returndds <- function(samples, meta_dir, tx2gene, dir){
  
  #List all possible comparisons available for experiment
  condition_pairs <- t(combn(levels(samples$condition), 2)) 
  
  #Print the options
  print("These are the available condition pairs")
  print(condition_pairs)
  
  #invisible(readline(prompt="Press [enter] to continue"))
  
  #Enter which comparison you would like to use 
  print("Please input which comparison you would like to choose (1, 2 or 3)")
  
  howzit <- scan(file = "", nmax = 1)
  
  #The particular comparison meta data (this cycles through the different ones)
  comparison <- samples[which(samples$condition == condition_pairs[howzit,]),]
  
  #Set rownmaes of comparison to sample ids
  rownames(comparison) <- comparison$sample
  
  #Get file paths for relevant comparison
  files_comparison <- file.path(dir, comparison$sample, "abundance.h5")
  
  #Set file names to be the sample ids
  names(files_comparison) <- comparison$sample
  
  #Make sure the comparison is a factor for Deseq2
  comparison$condition <- factor(comparison$condition)
  
  #Run tximport to get the kallisto files set txout to false to return gene level abundance
  txi.comparison <- tximport(files_comparison, type = "kallisto", txOut = FALSE, tx2gene = tx2gene)
  
  #Set colnames of the txi to be the sample ids
  colnames(txi.comparison$counts) <- rownames(comparison)
  
  #Create dds object for Deseq2
  dds <- DESeqDataSetFromTximport(txi.comparison, colData = comparison, design = ~condition)
  
  #Run DeSeq2
  dds <- DESeq(dds)
  
  
  #Can also sort by log2fold EMvNres <- subset(ddsresEMvN, padj < 0.05 & log2FoldChange > 0)
  print(praise("${Exclamation}! This is ${adjective}!"))
  #Return dds
  return(dds)
 
}
