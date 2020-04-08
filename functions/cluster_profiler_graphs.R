#####################################################################
'Library vibe'
#####################################################################

library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(forcats)
library(enrichplot)

rm(GO_enrichall)
#salmon <- GO_enrichall
#skyblue3 <- GO_enrichall
edo <- NULL
edo <- skyblue3



royalblue <- #GO_enrichall
  
  
#png(filename=names(normalised)[i], width = 700, height = 600, )
#dotplot(GO_enrichall, showCategory = 30)
#dev.off()
  enrichGO(gene = genes_GO, 			#Specify the gene set must be entrez ids
           universe = universe_entrez, #Universe refers to all the genes that it will be compared to. This would be a list of all the genes in Entrez IDs from either the annotation file, alternatively may be all the genes found on the probe.
           OrgDb = org.Hs.eg.db,	#The gene set must be annotated with the GO terms and requires a database to do so against. The correct database for the organsims which you're using is needed. Alex you need to download the org.Hs.eg.db database from bioconductor and use it instead
           ont = ont, 				#This specifies what ontology you wish to do enrichment on (BP= Biological process, CC= cellular component and MF= Molecular function)
           pAdjustMethod = "BH",	# this specifies the statistical method to perform on the genes. Benjamini and hochburg (BH) is prefered for GOEA (GO enrichment analysis)
           pvalueCutoff = 0.05,		# Cut off point for significant genes
           qvalueCutoff = 0.05, 
           minGSSize = 0, 
           maxGSSize = 2000,
           readable = TRUE)
  
edo@gene 
kegg <- enrichKEGG(gene = edo@gene, organism = 'hsa')
dotplot(edo)

enrichKEGG()

  
  
#####################################################################
'Prep Gene List File'
#####################################################################
#edo <- royalblue
#rm(GO_enrichall)
####################
#Get sig info from gene_info file
log2folds <- gene_info[gene_info$SYMBOL %in% edo@gene2Symbol,]
log2folds <- log2folds[,c(3,5)]

## assume that 1st column is ID
## 2nd column is fold change

## feature 1: numeric vector
geneList <- d[,2]

## feature 2: named vector
names(geneList) <- as.character(d[,1])

## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)

data(geneList, package="DOSE")
head(geneList)



#####################################################################
'Upset plot'
#####################################################################

upsetplot(edo, 15)
?upsetplot()

#####################################################################
'Bar Plot'
#####################################################################

library(enrichplot)
barplot(edo, showCategory=20)


#####################################################################
'Dot and Box'
#####################################################################
edo2 <- gseNCG(geneList, nPerm=10000)
dotplot(edo,       							#Drop True removes the significance values. Preferred when showing multiple plots
        showCategory = 20)
edo@ontology

update <- dropGO(edo, term = )

update$ term
edo@
  require(gridExtra)
?dropGO
require(clusterProfiler)


enrichMap(edo)
#####################################################################
'Gene-Concept Network'
#####################################################################
#bp2 <- simplify(edo, cutoff=0.7, by="p.adjust", select_fun=min)
simplify()


#edox
## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(edox, categorySize="pvalue", foldChange=geneList)
cnetplot(edox, foldChange=geneList, circular = T, colorEdge = TRUE, layout = "kk", showCategory = 10)
?cnetplot


require(DOSE)
require(enrichplot)
require(viridis)
data(geneList)
x = enrichDO(names(geneList)[1:100])
dotplot(x) + scale_color_viridis()







?cnetplot
dev.off()
par
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
dev.off()
grid.arrange(p1, p3, ncol=2)
test <- dropGO(GO_enrichall, term = 3)

?dropGO
# log2folds <- gene_info[gene_info$SYMBOL %in% GO_enrichall@gene2Symbol,]
# log2folds <- log2folds[,c(3,5)]
# log2folds <- log2folds[,c(2)]
# names(log2folds)[1] <- 'geneList'
# names(log2folds)[2] <- 'foldChange'


#####################################################################
'Heatmap'
#####################################################################

p1 <- heatplot(edox)
p2 <- heatplot(edox, foldChange=geneList)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
heatplot(edox, foldChange=geneList)


#####################################################################
'Enrichment Map'
######################################################################
edo$size <- edo$ Count

p1 <- emapplot(edox)
p2 <- emapplot(edo, pie_scale=1.5)
p3 <- emapplot(edo,layout="kk")
p4 <- emapplot(edo, pie_scale=1.5,layout="kk") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

library(clusterProfiler)
data(gcSample)
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
p1 <- emapplot(xx)
p2 <- emapplot(xx,legend_n=2) 
p3 <- emapplot(xx,pie="count")
p4 <- emapplot(xx,pie="count", pie_scale=1.5, layout="kk")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])


upsetplot(kk2) 

ridgeplot(edox)
upsetplot(edo)

#####################################################################
'Pubmed plot'
#####################################################################

terms <- edo$Description[1:5]
pmcplot(terms, 2010:2019, proportion=FALSE)


p <- pmcplot(terms, 2010:2017)

p2 <- pmcplot(terms, 2010:2019, proportion=FALSE)


plot_grid(p, p2, ncol=2)

p1 <- heatplot(edox)
p2 <- heatplot(edox, foldChange=geneList)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])









p1 <- emapplot(edo)
p2 <- emapplot(edo, pie_scale=1.5)
p3 <- emapplot(edo,layout="kk")
p4 <- emapplot(edo, pie_scale=1.5,layout="kk") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])






#####################################################################################################################################
#####################################################################################################################################
