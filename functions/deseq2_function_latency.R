#Deseq Latency Study

library(DESeq2)
library(tximport)



dir <- "~/doorknobdave/alecS/raw/latencystudy"
setwd(dir)
meta_dir <- "/home/studentsgh129/doorknobdave/alecS/raw"
tx2gene <- read.csv(file.path(meta_dir, "tx2gene94.txt"))

#Input meta data with accession numbers, cell type and time taken
samples <- read.csv(file.path(meta_dir, "meta.csv"), sep = ",", header = TRUE, stringsAsFactors = TRUE)


#Set condition for the run 
samples$condition <- samples$cell
samples$sample <- paste(samples$sample,samples$run_accession, sep = "_")

#Generate a list of all possible pairwise contrasts in my case this is only 3
condition_pairs <- t(combn(levels(samples$condition), 2))

for (howzit in 1:nrow(condition_pairs)){

#The particular comparison meta data (this cycles through the different ones)
comparison <- samples[which(samples$condition == condition_pairs[howzit,1] | samples$condition == condition_pairs[howzit,2]),]

#Set rownmaes of comparison to sample ids
rownames(comparison) <- comparison$sample

#Get file paths for relevant comparison
files_comparison <- file.path(dir, comparison$sample, "abundance.h5")

#Set file names to be the sample ids
names(files_comparison) <- comparison$sample

#Make sure the comparison is a factor for Deseq2
comparison$condition <- factor(comparison$condition)
tx2gene

file.exists(files_comparison)
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
#outdir <- "~/doorknobdave/alecS/results/latencyresults"

#Write to file
write.csv(sig_genes , paste(condition_pairs[howzit,1], condition_pairs[howzit,2], "significant.csv", sep = "_"))

}







#Main loop for all comparisions

howzit <- 2



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
  