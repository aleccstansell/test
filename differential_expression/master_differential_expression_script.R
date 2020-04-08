#########################################################################################################
#########################################################################################################
'Differential Expression Lists Master Script'
#########################################################################################################
#########################################################################################################

#########################################################################################################
'Preparation'
#########################################################################################################

#Load Packages
packages <- c("gplots", "ggplot2", "knitr", "limma", "reshape2", "dplyr", "RColorBrewer", "WGCNA", "praise", "Homo.sapiens","RColorBrewer")
lapply(packages, require, character.only = TRUE)

#Set Working directory 
out_dir <- "~/doorknobdave/alecS/alec_code/output"
meta_dir <- "~/doorknobdave/alecS/meta"
inputdir <- "~/doorknobdave/alecS/intermediary_data/count_matrices"

#Read meta data in 
meta_data <- read.table(file.path(meta_dir, "meta_for_coexpression_30.csv"), header = TRUE, sep = ",", col.names = c("sample_id", "condition", "patient"), row.names = NULL, colClasses = c("character", "factor", "factor"))

#Read in count data that is summarised to gene level see doc work flow for a list of script in order
raw_counts <- read.table(file.path(inputdir, "genelevelabundance_30.csv"), sep = ",", row.names=1)

#Running from 8 samples following outlier removal

#Load from coexpression output
setwd('/home/studentsgh129/doorknobdave/alecS/alec_code/coexpression/final_coexpression/data')
load("final_heat_normalised.RData")

#Set raw and meta data
raw_counts <- normalised_clean$input_counts
meta_data <- normalised_clean$`Input Metadata`


###################
#New
####################
#Load from coexpression output
setwd('/home/studentsgh129/doorknobdave/alecS/alec_code/coexpression/final_coexpression/results/Jan2020')
load("clean_all_samples2020.RData")
raw_countsnew <- clean_all_samples$input_counts
meta_datanew <- clean_all_samples$`Input Metadata`
samplesnew <- meta_datanew


#########################################################################################################
'Get Raw Counts Matrix'
#########################################################################################################

library(tximport)

##Set base directories 
work_dir <- "~/doorknobdave/alecS/raw/latencystudy/comparison"
meta_dir <-  "~/doorknobdave/alecS/raw/latencystudy/meta"
tx2gene <- read.csv(file.path("~/doorknobdave/alecS/meta", "tx2gene94.txt"))

#Read in sample meta data for file names
meta <- read.csv(file.path(meta_dir, "latency_meta.csv"), header = TRUE)

#Set rownmaes of comparison to sample ids
rownames(comparison) <- meta$sample_title

#Get file paths for relevant comparison
files_import <- file.path(work_dir, meta$sample, "abundance.h5")

#Name correctly
names(files_import) <- meta$sample_title

#Run tximport to get the kallisto files set txout to false to return gene level abundance
txi.object_healthy <- tximport(files_import, type = "kallisto", txOut = FALSE, tx2gene = tx2gene)

counts <- txi.object_healthy$counts
colnames(counts) <- meta$sample
 
#Write to file
write.csv(counts, file.path('~/doorknobdave/alecS/raw/latencystudy','Counts_latency.csv' ))

raw_counts <- counts
meta_data <- meta


#########################################################################################################
'LimmaR Differential Expression'
#########################################################################################################


############################################
#Load Functions
############################################

#Load Limma Script
source('~/doorknobdave/alecS/alec_code/coexpression/final_coexpression/limma_coexpression.R')

#Load Distribution Plot Script
source('~/doorknobdave/alecS/alec_code/coexpression/final_coexpression/distribution_plot.R')

############################################
#Run Limma using function
############################################
library(dplyr)

#Get gene lists for limmaR comparisons - in this case 3
CM_EM <- unlock_the_limma(raw_counts, meta_data)
CM_N <- unlock_the_limma(raw_counts, meta_data)
EM_N <- unlock_the_limma(raw_counts, meta_data)

#Generate List to store genes in 
limmasig_latency <- list(CM_EM, CM_N, EM_N)

#Name list items for easy storage
names(limmasig_latency) <- c("CM_EM", "CM_N", "EM_N")

############################################
#Save Files
############################################

#Set dir you want files to go
file_out <- "~/doorknobdave/alecS/alec_code/degs_output/latencystudy/limma"

#Save r object
save(limmasig_latency, file = file.path(file_out,"limmasig_genes_list.RData"))

#load(file.path(file_out,"limmasig_genes_list.RData"))

#Save all lists 
for(i in 1:3){
  
  name <- paste(names(limmasig_latency)[i], "limma_significant_no_dups.csv", sep = "_")
  
  write.csv(limmasig_latency[[i]], file.path(file_out, name))
  
  
}


#########################################################################################################
'DESeq2'
#########################################################################################################

##############################################
#Run DESeq2
##############################################

#Files needed for function of deseq which does tximport and gets list of all comparisons sig genes
#samples <-  read.table(file.path(meta_dir, "meta_for_coexpression_30.csv"), header = TRUE, sep = ",", col.names = c("sample_id", "condition", "patient"), row.names = NULL, colClasses = c("character", "factor", "factor"))
samples <- meta_data
#samples$sample <- samples$sample_title
tx2gene <- read.csv(file.path(meta_dir, "tx2gene94.txt"))
meta_dir <-  "~/doorknobdave/alecS/raw/latencystudy/meta"
dir <-"~/doorknobdave/alecS/raw/latencystudy/comparison"
file_out_deseq <- "~/doorknobdave/alecS/alec_code/degs_output/latencystudy/deseq2"


#Load function for deseq
scriptdir <- "~/alec/alec_code/functions"

#Load DE
source('~/doorknobdave/alecS/alec_code/functions/deseq_only_no_loop.R')
#source('~/doorknobdave/alecS/alec_code/functions/deseq_return_dds.R')

#Run Deseq2 returning alist of sig genes and a dds object
EMvsCM_sig_genes <- deseqonly(samples = samples, meta_dir = meta_dir, tx2gene = tx2gene, dir = dir)
NvsCM_sig_genes <- deseqonly(samples = samples, meta_dir = meta_dir, tx2gene = tx2gene, dir = dir)
NvsEM_sig_genes <- deseqonly(samples = samples, meta_dir = meta_dir, tx2gene = tx2gene, dir = dir)

#Select dds object for saving purposes
ddsEMvsCM <- EMvsCM_sig_genes[[2]]
ddsNvsCM <- NvsCM_sig_genes[[2]]
ddsNvsEM <- NvsEM_sig_genes[[2]]

#Save DDS Objects
deseq_ddsobjects <- list(EMvsCM_sig_genes, NvsCM_sig_genes, NvsEM_sig_genes)
names(deseq_ddsobjects) <- c("EMvsCM_sig_genes", "NvsCM_sig_genes", "NvsEM_sig_genes")
save(deseq_ddsobjects, file = file.path(file_out_deseq,"DDSobjects_deseq.RData"))


#######################################
#New 8
#####################################
EMvsCM_sig_genes_new <- deseqonly(samples = samplesnew, meta_dir = meta_dir, tx2gene = tx2gene, dir = dir)
dds_EMvsCMnew <-  EMvsCM_sig_genes_new[[2]]

NvsCM_sig_genes_new <- deseqonly(samples = samplesnew, meta_dir = meta_dir, tx2gene = tx2gene, dir = dir)
dds_NvsCMnew <-  NvsCM_sig_genes_new[[2]]

NvsEM_sig_genes_new <- deseqonly(samples = samplesnew, meta_dir = meta_dir, tx2gene = tx2gene, dir = dir)
dds_NvsEMnew <- NvsEM_sig_genes_new[[2]]



#Save DDS Objects
file_out_deseq <- "~/doorknobdave/alecS/alec_code/degs_output/deseq2_degs"
final_8_ddsobjects <- list(dds_EMvsCMnew, dds_NvsCMnew, dds_NvsEMew)
names(final_8_ddsobjects) <- c("dds_EMvsCMnew", "dds_NvsCMnew", "dds_NvsEnew")
save(final_8_ddsobjects, file = file.path(file_out_deseq,"final_8_ddsobjects.RData"))



#Return dds
# ddsEMvsCM <- deseq_returndds(samples = samples, meta_dir = meta_dir, tx2gene = tx2gene, dir = dir)
# ddsNvsCM <- deseq_returndds(samples = samples, meta_dir = meta_dir, tx2gene = tx2gene, dir = dir)
# ddsNvsEM <- deseq_returndds(samples = samples, meta_dir = meta_dir, tx2gene = tx2gene, dir = dir)

#Get results
ddsres <- results(ddsEMvsCM)
nrow(ddsres)
#Order by dds p adjusted
sig_genes <- subset(ddsres, padj < 0.05)
nrow(sig_genes)


##############################################
#Save Deseq2 to file or load it
##############################################

#List of deseq
deseq_gene_list <- list(EMvsCM, NvsCM, NvsEM)
names(deseq_gene_list) <- c("EMvsCM", "NvsCM", "NvsEM")

#Set dir you want files to go
file_out_deseq <- "~/doorknobdave/alecS/alec_code/degs_output/deseq2_degs"

#Save Significantly Expressed Gene Object
save(deseq_gene_list, file = file.path(file_out_deseq,"deseq_gene_list.RData"))

#Save DDS Objects
deseq_ddsobjects <- list(ddsEMvsCM, ddsNvsCM, ddsNvsEM)
names(deseq_gene_list) <- c("ddsEMvsCM", "ddsNvsCM", "ddsNvsEM")
save(deseq_gene_list, file = file.path(file_out_deseq,"DDSobjects_deseq.RData"))

#Load Deseq2 file
load(file.path(file_out_deseq,"deseq_gene_list.RData"))
deseq_gene_list$NvsEM

############################################
#Save Files
############################################


deseqfinal8signficant <- list(EMvsCM_sig_genes_new[[1]], NvsCM_sig_genes_new[[1]], NvsEM_sig_genes_new[[1]])
names(deseqfinal8signficant) <- c("EMvsCM", "NvsCM", "NvsEM")
save(deseqfinal8signficant, file = file.path(file_out_deseq,"final8significantgenes.RData"))

#############################
#Healthy
###########################3

deseqhealthy <- list(EMvsCM_sig_genes[[1]], NvsCM_sig_genes[[1]], NvsEM_sig_genes[[1]])
names(deseqhealthy) <- c("EMvsCM", "NvsCM", "NvsEM")
save(deseqhealthy, file = file.path(file_out_deseq,"sig_latency_genes.RData"))
load(file = file.path(file_out_deseq,"sig_healthy_cells.RData"))


#Save all lists 
for(i in 1:3){
  
  name <- paste(names(deseqhealthy)[i], "deseq2_significant_no_dups.csv", sep = "_")
  
  write.csv(deseqhealthy[[i]], file.path(file_out_deseq, name))
}  
  

##############################################
#Select Upregulated Down regulated etc
##############################################

sig <- deseq_gene_list$EMvsCM
sig <- deseq_gene_list$NvsCM
sig <- deseq_gene_list$NvsEM

#Choose whether you want upregulated downregulated or all differentially expressed genes
sig <- sig[,-1]
colnames(sig) <- c("gene_names", "log2fold")


#Print choices
print("Choose whether you would like to use downregulated, upregulated or all DE genes")
invisible(readline(prompt="Press [enter] to continue"))
print("Type up for upregulated, down for downregulated and all for all Differentially expressed genes")
reg <- scan(file = "", what = "character", nmax = 1)

reg <- 'up'
i <- 1


#######################################################################################################################
'Count Number of Genes upregulated and downregulated'
#######################################################################################################################

#####################################################
#Count DESeq2
#####################################################

for(i in 1:3){
  
  #Upregulated
  df <- deseqhealthy[[i]]
  df <- subset(df, log2FoldChange > 0)
  print(paste(names(deseqfinal8signficant[i]), nrow(df), "Upregulated genes selected", sep = " "))
  
  #Down
  df <- deseqhealthy[[i]]
  df <- subset(df, log2FoldChange < 0)
  print(paste(names(deseqfinal8signficant[i]), nrow(df), "Downregulated genes selected", sep = " "))
  
  
  
}

#####################################################
#Count Sleuth
#####################################################

nrow(list_sleuth$EMvsCM_sig_sleuth)
nrow(list_sleuth$CMvsN_sig_sleuth)
nrow(list_sleuth$EMvsN_sig_sleuth)

#####################################################
#Count Limma
#####################################################
for(i in 1:3){
  
  #Upregulated
  df <- limmasig_genes[[i]]
  df <- subset(df, logFC > 0)
  print(paste(names(limmasig_genes[i]), nrow(df), "Upregulated genes selected", sep = " "))
  
  #Down
  df <- limmasig_genes[[i]]
  df <- subset(df, logFC < 0)
  print(paste(names(limmasig_genes[i]), nrow(df), "Downregulated genes selected", sep = " "))
  
  
  
}

# 
# dds <- final_8_ddsobjects
# results(dds)
# 
# 
# nrow(df)
# 
# 
# if(reg == "up") {
#   df <- subset(sig, log2FoldChange > 1)
#   df <- df[order(df$log2FoldChange, decreasing = FALSE),]
#   print("Upregulated genes selected")
# } else if(reg == "down"){
#   df <- subset(sig, log2FoldChange < 1)
#   df <- df[order(df$log2FoldChange, decreasing = TRUE),]
# 
# 
# 
# 
# #If statement to select specific genes needed 
# if(reg == "up") {
#   df <- subset(sig, log2FoldChange > 1)
#   df <- df[order(df$log2FoldChange, decreasing = FALSE),]
#   print("Upregulated genes selected")
# } else if(reg == "down"){
#   df <- subset(sig, log2FoldChange < 1)
#   df <- df[order(df$log2FoldChange, decreasing = TRUE),]
#   print("Downregulated genes selected")
# } else if(reg == "all"){
#   df <- sig 
#   df <- df[order(df$log2FoldChange, decreasing = TRUE),]
#   print("All genes selected")
#   
# }
# nrow(df)
# nrow(sig)

#######################################################################################################################
'Sleuth'
#######################################################################################################################

#Go to this script and run - this is not a function

#Use commander script for analysis
source('~/doorknobdave/alecS/alec_code/sleuth/commander_master_sleuth_30.R')



##############################################
#Venn Diagram
##############################################

#Change wd to output venn diagrams
setwd("~/doorknobdave/alecS/alec_code/degs_output/venn_diagrams")

#Load Venn Script
source('~/doorknobdave/alecS/alec_code/functions/venn_diagram.R')

#Load Deseq gene list
file_out_deseq <- "~/doorknobdave/alecS/alec_code/degs_output/HIV_control/deseq2_degs"
load(file.path(file_out_deseq,"final8significantgenes.RData"))

#Load Limma gene lists
file_out_limma <- "~/doorknobdave/alecS/alec_code/degs_output/limma_degs"
load(file.path(file_out_limma,"limmasig_genes_list.RData"))

#Load sleuth gene lists
file_out_sleuth <- "~/doorknobdave/alecS/alec_code/degs_output/sleuth_degs"
load(file.path(file_out_sleuth,"sleuthsig_genes_list.RData"))

############################################################################################
#Generate Venn Diagrams using venn script
############################################################################################

getwd()
setwd('~/doorknobdave/alecS/alec_code/degs_output/venn_diagrams')

#Venn for CM EM
vennit(list_sleuth$EMvsCM_sig_sleuth$target_id, rownames(deseqfinal8signficant$EMvsCM), limmasig_genes$CM_EM$Gene_Name, "sleuth", "DESeq2", "limma", "EMvsCM_venn.png")

#Venn for CM N
vennit(list_sleuth$CMvsN_sig_sleuth$target_id, rownames(deseqfinal8signficant$NvsCM), limmasig_genes$CM_N$Gene_Name, "sleuth", "DESeq2", "limma", "CMvsN_venn.png")

#Venn for EM N
vennit(list_sleuth$EMvsN_sig_sleuth$target_id, rownames(deseqfinal8signficant$NvsEM), limmasig_genes$EM_N$Gene_Name, "sleuth", "DESeq2", "limma", "EMvsN_venn.png")


list_sleuth$EMvsCM_sig_sleuth

#Check all genes are signficant Benjamini Hochburg adjusted less than .05
all(list_sleuth$EMvsN_sig_sleuth$qval < 0.05)
all(list_sleuth$EMvsN_sig_sleuth$qval < 0.05)
all(list_sleuth$EMvsN_sig_sleuth$qval < 0.05)


sleuth_all <- cbind(list_sleuth$EMvsCM_sig_sleuth$target_id, list_sleuth$CMvsN_sig_sleuth$target_id, list_sleuth$EMvsN_sig_sleuth$target_id)

############################################################################################
"Get Genes That Are Differentially Expressed In ALl 3 Methods"
############################################################################################

#EMvsCM All Sig For ALl Methods
final <- NULL
deseq_sleuth <- NULL
all <- NULL

#Get lists and accumulate
rownames(list_sleuth$EMvsCM_sig_sleuth) <- list_sleuth$EMvsCM_sig_sleuth$target_id
rownames(limmasig_genes$CM_EM)
deseq_sleuth <- merge(as.data.frame(deseqfinal8signficant$EMvsCM), as.data.frame(list_sleuth$EMvsCM_sig_sleuth), by="row.names", sort=FALSE)
rownames(deseq_sleuth) <- deseq_sleuth$Row.names
all <- merge(deseq_sleuth, as.data.frame(limmasig_genes$CM_EM),  by="row.names", sort=FALSE)
rownames(all) <- all$Row.names
all$Row.names <- NULL
final <- merge(all, as.data.frame(normalised_clean$variance_stabilised_counts),  by="row.names", sort=FALSE)
final$logFC <- -(final$logFC)
head(final)
EMvsCM_all <- final


#CMvsN All Sig For ALl Methods
final <- NULL
deseq_sleuth <- NULL
all <- NULL

#Get lists and accumulate
rownames(list_sleuth$CMvsN_sig_sleuth) <- list_sleuth$CMvsN_sig_sleuth$target_id
rownames(limmasig_genes$CM_N)
deseq_sleuth <- merge(as.data.frame(deseqfinal8signficant$NvsCM), as.data.frame(list_sleuth$CMvsN_sig_sleuth), by="row.names", sort=FALSE)
rownames(deseq_sleuth) <- deseq_sleuth$Row.names
all <- merge(deseq_sleuth, as.data.frame(limmasig_genes$CM_N),  by="row.names", sort=FALSE)
rownames(all) <- all$Row.names
all$Row.names <- NULL
final <- merge(all, as.data.frame(normalised_clean$variance_stabilised_counts),  by="row.names", sort=FALSE)
final$logFC <- -(final$logFC)
NvsCM_all <- final



#EMvsCM All Sig For ALl Methods
final <- NULL
deseq_sleuth <- NULL
all <- NULL

#Get lists and accumulate
rownames(list_sleuth$EMvsN_sig_sleuth) <- list_sleuth$EMvsN_sig_sleuth$target_id

deseq_sleuth <- merge(as.data.frame(deseqfinal8signficant$NvsEM), as.data.frame(list_sleuth$EMvsN_sig_sleuth), by="row.names", sort=FALSE)
rownames(deseq_sleuth) <- deseq_sleuth$Row.names
all <- merge(deseq_sleuth, as.data.frame(limmasig_genes$EM_N),  by="row.names", sort=FALSE)
rownames(all) <- all$Row.names
all$Row.names <- NULL
final <- merge(all, as.data.frame(normalised_clean$variance_stabilised_counts),  by="row.names", sort=FALSE)
final$logFC <- -(final$logFC)
head(final)
NvsEM_all <- final

#################Save results#####################

results_all_methods <- list(EMvsCM_all, NvsCM_all, NvsEM_all)
names(results_all_methods) <- c("EMvsCM_all", "NvsCM_all", "NvsEM_all")
setwd("~/doorknobdave/alecS/alec_code/degs_output/occuring_all")
save(results_all_methods, file = "list_genes_occuring_in_all_methods.RData")
load("list_genes_occuring_in_all_methods.RData")

file_out <- "~/doorknobdave/alecS/alec_code/degs_output/occuring_all"
#Save all lists 
for(i in 1:3){
  
  name <- paste(names(results_all_methods)[i], "all_3_methods.csv", sep = "_")
  
  write.csv(results_all_methods[[i]], file.path(file_out, name))
  
  
}






list_sleuth$EMvsCM_sig_sleuth

deseq_gene_list$EMvsCM
  
limmasig_genes$CM_EM
  


