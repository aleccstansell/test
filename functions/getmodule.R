#####################################################################################################################################
#Extract Genes from modules function
#####################################################################################################################################

'Input - list of genes from cytoscape with attached column module containing module number'
'Output - Genes in that module'

getmodule <- function(gene_list, modulenum){
  gene2 <- gene_list[which(gene_list$module == modulenum),]
  #gene1 <- gene_list$name
  return(gene2)
}





####################################

# setwd("/home/rembrandt/alec/results")
# input <- read.csv("SOM_cluster_400_EM_CM_deseq.csv")
# 
# head(input)
# 
# input$module <- input$cluster_SOM 
# 
# mod3 <- getmodule(input, modulenum = 2)
# 
# input$cluster_SOM <- as.factor(input$cluster_SOM)
# 
# input$cluster_kmeans <- as.factor(input$cluster_kmeans)
# 

