#Almagamate Module GO terms
#Joins together a big table of all GO terms
##Set base directories 

#file_location <- "~/doorknobdave/alecS/alec_code/coexpression/final_coexpression/results/30157_genes/cluster_profiler/NCMwgcna"

#out_dir <- "~/doorknobdave/alecS/alec_code/coexpression/final_coexpression/results/30157_genes/cluster_profiler"

#name <- "EMvsCM"



almalgamate_GO <- function(file_location, out_dir){
  
  setwd(file_location)
  
  #Create object file names containing file paths for all the files 
  file.names <- file.path(file_location, list.files(file_location))
  
  #initialise
  master <- NULL
  #Loop through files to add all module GO enriched data together
  for (i in 1:length(file.names)){
    #Read file particular file in
    GO <- read.table(file.names[i],header=TRUE, sep=",", stringsAsFactors=FALSE) 
  
      if(nrow(GO) > 0){
      #Get list of all module names and add them to each row
      GO$module <- noquote(unlist(strsplit(list.files(file_location), split = "_GOenrichallprocess.csv", fixed = TRUE))[i])
      master <- rbind(master, GO)
      
    } else {
      print(paste(noquote(unlist(strsplit(list.files(file_location), split = "_GOenrichallprocess.csv", fixed = TRUE))[i]), "has no data", sep = " "))
    }
    
    
  }
 
  #Write to file
  write.csv(master, file.path(out_dir, paste(name, 'WGCNA_all_GO_TERMS.csv', sep = "")))
  
  
  
}

# work_dir <- "~/doorknobdave/alecS/alec_code/coexpression/final_coexpression/results/30157_genes/cluster_profiler/all_methods/NvsEM_wgcna_allmethods/modules"
# 
# getwd()
# 
# #work_dir <- "~/doorknobdave/alecS/alec_code/clustering/wgcna/results/cluster_profiler/NEMwgcna"
# setwd(work_dir)
# #Create object file names containing file paths for all the files
# file.names <- file.path(work_dir, list.files(work_dir))
# 
# #initialise
# master <- NULL
# #Loop through files to add all module GO enriched data together
# for (i in 1:length(file.names)){
#   #Read file particular file in
#   GO <- read.table(file.names[i],header=TRUE, sep=",", stringsAsFactors=FALSE)
# 
# 
#   if(nrow(GO) > 0){
#         #Get list of all module names and add them to each row
#     GO$module <- noquote(unlist(strsplit(list.files(work_dir), split = "_GOenrichallprocess.csv", fixed = TRUE))[i])
#     master <- rbind(master, GO)
# 
#   } else {
#     print(paste(noquote(unlist(strsplit(list.files(work_dir), split = "_GOenrichallprocess.csv", fixed = TRUE))[i]), "has no data", sep = " "))
#   }
# 
# 
# }
# 
# out_dir <- "~/doorknobdave/alecS/alec_code/coexpression/final_coexpression/results/30157_genes/cluster_profiler/all_methods/"
# 
# 
# write.csv(master, file.path(out_dir, 'N_EM_WGCNA_all_GO_TERMS.csv'))


#   source(file.path(scriptdir, "convert_between_ids.R"), echo = TRUE)
# 
# 
# genes <- convert_between_ids(getmodule(gene_list,50),input_convert)
# 
# genes <- getmodule(gene_list, 50)
# genes <- getentrez(genes)
# 
# clusterprop <- clusterit(genes, universe_entrez, 'all')
# 
# getwd()
# 
# GO_enrichall <- enrichGO(gene = genes, 			#Specify the gene set must be entrez ids
#                          universe = universe_entrez, #Universe refers to all the genes that it will be compared to. This would be a list of all the genes in Entrez IDs from either the annotation file, alternatively may be all the genes found on the probe.
#                          OrgDb = org.Hs.eg.db,	#The gene set must be annotated with the GO terms and requires a database to do so against. The correct database for the organsims which you're using is needed. Alex you need to download the org.Hs.eg.db database from bioconductor and use it instead
#                          ont = "all", 				#This specifies what ontology you wish to do enrichment on (BP= Biological process, CC= cellular component and MF= Molecular function)
#                          pAdjustMethod = "BH",	# this specifies the statistical method to perform on the genes. Benjamini and hochburg (BH) is prefered for GOEA (GO enrichment analysis)
#                          pvalueCutoff = 0.05,		# Cut off point for significant genes
#                          qvalueCutoff = 0.05, 
#                          minGSSize = 0, 
#                          maxGSSize = 2000,
#                          readable = TRUE)		
# 
# 
# 
# GO_enrichall$ Description
# 
# #By setting reable true will attempt to write the entrez IDs back into gene symbols. Was placed as false as many genes aren't annotated with a gene symbol yet.
# write.table(GO_enrichall, file = paste("module_", module, "_GOenrichallprocess.csv", sep = ""), sep = ",")
# 
# 
# 
# input_convert <- c('ensembl_gene_id', 'hgnc_symbol')
# getmodule(genes_list = gene_list, modulenum = 50)
# 
# 
# write.csv(master, 'all_modules_GO_terms.csv')
# 
# getwd()
# 
# #Check number of duplicated values
# 
# table(duplicated(master$Description))
# 
# 
# 
# master
# 
# all(unique(master$module))
# 
# table(master$Description)
# 
# 
# 
# library(clusterProfiler)
# 
# 
# 
# 
# 
# 














# 
# 
# 
# table(duplicated(master$Description))
# 
# 
# duplicated_values <- master[which(duplicated(master$Description)),]
# 
# length(duplicated_values$Description)
# 
# 
# 
# write.csv(duplicated_values, "duplicate_GO_terms.csv")
# 
# 
# 
# sprintf(paste(as.character(c))
# 
# a <- length(unique(master$Description))
# b <- length(master$Description)
# 
# c <- b-a
# 
# 
# master$Description
# 
# master$module
# 
# 
# 
# 
# master
# unique(GO)
# 
# length(master$ID) == length(unique(master$ID))
# 
# table(master$module)
# 
# 
# tail(master)
# 
# 
# 
# master$module
# 
# #Write to csv
# write.csv(master, "all_GO_for_deseq_coexpression.csv")
