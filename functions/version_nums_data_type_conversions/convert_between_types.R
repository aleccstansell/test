#Install Packages
#biomart

#Load Packages

library(biomaRt)


#Load in data to be converted

inputdir <- "~/doorknobdave/alecS/intermediary_data/Kallisto_all"
raw_counts <- read.csv(file.path(inputdir, "kallisto_results_68.csv"), row.names=1)
head(raw_counts)
dim(raw_counts)

#Create object gene_list with the ensembl gene ids
gene_list <-rownames(raw_counts) 

#Download biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

#View different types of attributes that can be interchanged
attributes(mart)

#gene_list=unique(sort(gather(df[-1],"","genes")$genes))

# Test using smaller gene list
genes <-gene_list[1:50]

#Get Biomart gene ids using attributes = which ids you want to extract, filters only entrez ids that have ensembl identifiers, values - the list of genes that you want to get the attributes for, bmHeader - whether there is a header and the object to use - mart)

#Example used to get entrez and ensembl for cluster profiler
#entrez=getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = "ensembl_transcript_id", values = genes, bmHeader = F, mart = mart)

entrez=getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = "ensembl_transcript_id", values = gene_list, bmHeader = F, mart = mart)



#Change colnames to reflect entrez gene ids

colnames(entrez)=c("ensembl_gene_id","entrezgene")

#List of entrez gene ids for the input list of ensembl genes is give by the following column
genes_for_GO <- entrez$entrezgene

genes_for_GO <- na.omit(genes_for_GO)

#Get Biomart Universe genes to compare all GO terms against
universe_entrez = getBM(attribute = c("entrezgene"), mart=mart)

universe_entrez <- universe_entrez$entrezgene
#Ensure correct class type for lists (ie not intergers - must be characters)

genes_for_GO <- as.character(genes_for_GO)


universe_entrez <- as.character(universe_entrez)

class(universe_entrez)

















