library(praise)
library(tximport)
library(DESeq2)

#Initialise
sig <- NULL

deseqit <- function(samples, meta_dir, tx2gene, dir){
  
  #Generate a list of all possible pairwise contrasts in my case this is only 3
  condition_pairs <- t(combn(levels(samples$condition), 2))          
  
  #Main loop for all comparisions
  
  for (howzit in 1:nrow(condition_pairs)){
    
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
    
    #Get results
    ddsres <- results(dds)
    
    #Order by dds p adjusted
    sig_genes <- subset(ddsres, padj < 0.05)
    
    #Can also sort by log2fold EMvNres <- subset(ddsresEMvN, padj < 0.05 & log2FoldChange > 0)
    print(praise("${Exclamation}! This is ${adjective}!"))
    
    #Sig genes
    sig <- union(rownames(sig_genes), sig)
    
    print(praise("${Exclamation}! This is ${adjective}!"))
    
  }
  
  return(sig)
  
}







