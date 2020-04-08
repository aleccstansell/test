#####################################################################################################################################
#Get entrez function - just input the genes and it will output entrez values
#####################################################################################################################################

'GET ENTREZ FUNCTION'
'INPUT - Ensembl IDS to be converted to entrez'
'OUTPUT - genes_for_GO as entrez IDS '

#Install Biomart
#BiocManager::install("biomaRt")

#Load biomart
library(biomaRt)

#Load mart object with the data
#mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
mart = useDataset("hsapiens_gene_ensembl",mart=ensembl)
# getwd()
# save(mart, list = character(0), file = "martobject.RData")
# load(martobject.Rdata)
# getwd()
# test <- read.table(file.path(outdir, "martobject.RData"))
# load(martobject.RData, envir = parent.frame(), verbose = FALSE)
# setwd("/home/studentsgh129/doorknobdave/alecS/intermediary_data/cytoscape_out/clusterout")
# load("martobject.RData")


genes <- genes_module$ENSEMBL
###############################################################################
#Get entrez function
##############################################################################
getentrez <- function(genes){
  
  #Get Biomart gene ids using attributes = which ids you want to extract, filters only entrez ids that have ensembl identifiers, values - the list of genes that you want to get the attributes for, bmHeader - whether there is a header and the object to use - mart)
  entrez=getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters = "ensembl_gene_id", values = genes, bmHeader = T, mart = mart)
  
  
  #Change colnames to reflect entrez gene ids
  colnames(entrez)=c("ensembl_gene_id","entrezgene_id")
  
  #List of entrez gene ids for the input list of ensembl genes is give by the following column
  genes_for_GO <- entrez$entrezgene
  
  #Remove na values
  genes_for_GO <- na.omit(genes_for_GO)
  
  #Ensure correct class type for lists (ie not intergers - must be characters)
  genes_for_GO <- as.character(genes_for_GO)
  
}



#####################################################################################################################################
#####################################################################################################################################
