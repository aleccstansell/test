#################################################################################################################################
#Convert between ids
#################################################################################################################################

#Load library
library(biomaRt)

#Create biomart object from ensembl
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

"For list of attributes to convert between see attributes(mart)"

#   input_type <- 'ensembl_gene_id'
#   input_genes <- genesin
#   type_to_convert_to <- 'hgnc_symbol'
#   input_convert <- c('hgnc_symbol','ensembl_gene_id')



convert_between_ids <- function(input_genes, input_convert){
  #Get Biomart gene ids using attributes = which ids you want to extract, filters only entrez ids that have ensembl identifiers, values - the list of genes that you want to get the attributes for, bmHeader - whether there is a header and the object to use - mart)
  ids=getBM(attributes = input_convert, filters = input_convert[1], values = input_genes, bmHeader = T, mart = mart)
  
  #Change colnames to reflect entrez gene ids
  colnames(ids)=input_convert
  
  output_genes <- as.character(ids[,2])
  return(output_genes)
  
}
# 
# 
# convert_between_ids(interest, input_convert)
# 
# output_genes
# 
# 
# 
# inputdir <- "~/doorknobdave/alecS/results/lists_deseq"
# 
# file <- "CMvsNpadj0.05.csv"
# 
# for (file in list.files(inputdir)){
#   
#   comp <- read.csv(file.path(inputdir, file))
#   comp$HGNC <- convert_between_ids(comp$X, input_convert)
#   write.csv(comp, file)
#   
#   
# }
# 
# 
# genes <- read.csv(file.path(inputdir, "EMvCMpadj0.05"))
# 
# files
# 
# 
# 
# 
# mod50genes <- convert_between_ids(input_genes = mod50, input_convert = input_convert)
# 
# output
# inputdir <- "~/doorknobdave/alecS/intermediary_data"
# 
# genes <- read.csv(file.path(inputdir, "EMvCMpadj0.05"))
# 
# genes$HGCN <- convert_between_ids(genes$names, input_convert)
# getwd()
# 
# write.csv(genes, 'hgcn_EMvCMpadj0.05.csv')
# 
# 
# 
# 
# convert_between_ids()
# 
#################################################################################################################################
#################################################################################################################################