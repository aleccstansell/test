#####################################################################################################################################
#Cluster Profiler Function
#####################################################################################################################################

'Inputs for cluster it are: genes_GO = entrez gene ids, universe_entrez - to be comapred and ont to be used BP, CC, ALL or MF'

library(clusterProfiler)

# genes_GO = genes_entrez
# universe_entrez = universe_entrez
# ont = 'all'
# 
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

#Get Biomart Universe genes to compare all GO terms against
universe_entrez = getBM(attribute = "entrezgene_id", mart=mart)
universe_entrez <- universe_entrez$entrezgene_id
universe_entrez <- as.character(universe_entrez)





clusterit <- function(genes_GO, universe_entrez, ont){
  #Run GO enrichmentfor Biological Process
  GO_enrichall <- enrichGO(gene = genes_GO, 			#Specify the gene set must be entrez ids
                           universe = universe_entrez, #Universe refers to all the genes that it will be compared to. This would be a list of all the genes in Entrez IDs from either the annotation file, alternatively may be all the genes found on the probe.
                           OrgDb = org.Hs.eg.db,	#The gene set must be annotated with the GO terms and requires a database to do so against. The correct database for the organsims which you're using is needed. Alex you need to download the org.Hs.eg.db database from bioconductor and use it instead
                           ont = ont, 				#This specifies what ontology you wish to do enrichment on (BP= Biological process, CC= cellular component and MF= Molecular function)
                           pAdjustMethod = "BH",	# this specifies the statistical method to perform on the genes. Benjamini and hochburg (BH) is prefered for GOEA (GO enrichment analysis)
                           pvalueCutoff = 0.05,		# Cut off point for significant genes
                           qvalueCutoff = 0.05, 
                           minGSSize = 0, 
                           maxGSSize = 2000,
                           readable = TRUE)		#By setting reable true will attempt to write the entrez IDs back into gene symbols. Was placed as false as many genes aren't annotated with a gene symbol yet.
  
  
  write.csv(data.frame(GO_enrichall), file.path("cluster_profiler", "module", paste(module, "_all_GO_terms", ".csv", sep = "")))
  
  save(GO_enrichall, file = file.path("cluster_profiler", "object", paste(module, "_cluster_profiler_object", ".RData", sep = "")))
 

  
  
}


# 
# png(filename=names(normalised)[i], width = 700, height = 600, )
# sampleHeatmap(normalised[[i]], names(normalised)[i])
# dev.off()
# 
# list_GO[[]]
# 
# png(filename=names(normalised)[i], width = 700, height = 600, )
# dotplot(GO_enrichall, showCategory = 30)
# dev.off()
# 
# invisible(readline(prompt="Press [enter] to continue"))
# 
# library(DOSE)
# 
# compareCluster()
# ?compareCluster()
# 
# enrichmap(GO_enrichall, 							
#           drop=FALSE, 								#Drop True removes the significance values. Preferred when showing multiple plots
#           showCategory = 12)							#Shows an interconnected map 
# 
# 
# dotplot(list_GO[[1]][[1]][[1]][[1]][[1]][[2]], showCategory = 30)
# 
# 
# 
# barplot(GO_enrichall, 							
#         drop=FALSE, 								#Drop True removes the significance values. Preferred when showing multiple plots
#         showCategory = 20)	
# 


#####################################################################################################################################
#####################################################################################################################################
