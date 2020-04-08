# Code to change ensembl gene ids to entrez ids

#Load Packages

library(biomaRt)
library(biomartr)
library(tidyr)
library(forcats)

install.packages("biomartr")
#Set working directory

work_dir <- "~/doorknobdave/final/results/deseqmarch2019/CMvsN"
setwd(work_dir)


#Read in dataset from Deseq2 containing ensembl gene ids

df <- read.csv("CMvsNpadj0.05.csv", sep=",", header = T, stringsAsFactors = F, strip.white = T)
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)

#Load mart object from biomart

#mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

mart = useDataset("hsapiens_gene_ensembl",mart=ensembl)

#Create object gene_list with the ensembl gene ids
gene_list=df$names

# gene_list=unique(sort(gather(df[-1],"","genes")$genes))

#Get Biomart gene ids using attributes = which ids you want to extract, filters only entrez ids that have ensembl identifiers, values - the list of genes that you want to get the attributes for, bmHeader - whether there is a header and the object to use - mart)

entrez=getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = "ensembl_gene_id", values = gene_list, bmHeader = T, mart = mart)

#Change colnames to reflect entrez gene ids

colnames(entrez)=c("ensembl_gene_id","entrezgene")

#List of entrez gene ids for the input list of ensembl genes is give by the following column
entrez$entrezgene

#See which IDs the bitr function can convert to and from
keytypes(org.Hs.eg.db)

#To get all BM entrez ids for use in clusterprofiler
all_entrez = getBM(attribute = c("entrezgene"), mart=mart)



entrez$entrezgene
colnames(test)=c("ensembl_gene_id","entrezgene")


test$ncbi_gene_id=as.character(test$ncbi_gene_id)

new_df=df
new_df[-1] <- lapply(new_df[-1], function(x) lvls_revalue(factor(x, levels = test$hgnc_symbol), test$ncbi_gene_id))