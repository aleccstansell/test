
###############################
#Install Packages
BiocManager::install("sva")

#Load packages
library(sva)

#Set dds
dds <- deseq_ddsobjects$NvsEM_sig_genes[2][[1]]
res_og <- results(dds)






file_out <- "~/doorknobdave/alecS/alec_code/degs_output/healthy_cells/deseq2"
setwd(file_out)
#Create coldata
coldata <- meta_data$sample
coldata <- meta_data$sample
coldata <- data.frame(coldata)
dds <- DESeqDataSetFromTximport(txi.comparison, colData = comparison, design = ~condition)


EMvN_raw <- raw_counts[, c(1, 3, 4,6, 7, 9)]
colnames(EMvN)
EMvN$run_accession <- NULL
rownames(EMvN) <- EMvN$sample


EMvN <- meta_data[ which(meta_data$cell=='EM'|meta_data$cell=='N'), ]
EMvN$run_accession <- NULL
rownames(EMvN) <- EMvN$sample



colnames(EMvN_raw) <- NULL

for(i in 1:6){
  EMvN_raw[,i] <- as.integer(EMvN_raw[,i])
  
}


rownames(comparison)
colnames(EMvN_raw)
as.integer()
EMvN 
# dds <-
ncol(EMvN_raw)
DESeq(dds)
dds <- DESeqDataSetFromMatrix(EMvN_raw, colData = EMvN, design = ~condition)

ncol(EMvN_raw)
dat <- counts(dds)
counts(dds)
#Data 
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv=2)
svseq$sv

getwd()
png("sva_plot", 1200, 1000, pointsize=5)
par(mfrow=c(2,1),mar=c(3,5,3,1))
stripchart(svseq$sv[,1] ~ dds$condition,vertical=TRUE,main="SV1")
abline(h=0)
stripchart(svseq$sv[,2] ~ dds$condition,vertical=TRUE,main="SV2")
abline(h=0)
dev.off()


ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + condition
ddssva <- DESeq(ddssva)
ddssva_res <- results(ddssva)
sig_genes <- subset(ddssva_res, padj < 0.05)
nrow(sig_genes)


rld <- rlogTransformation(ddssva)
pcaData <- plotPCA(rld, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

library(ggplot2)

png("adjusted.png", 1200, 1000, pointsize=5)
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size =10) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
dev.off()
getwd()



countmatrix<-read.csv(file="C:\\Kevin\\UniqueCounts.csv" ,header=TRUE, row.names=1)
as.matrix(countmatrix[,-1])
colnames(countmatrix) <- NULL
dds=DESeqDataSetFromMatrix(countData=countmatrix, colData=coldata, design= ~condition)
dds$test <-factor(paste0(dds$species,dds$treatment))
design(dds) <- ~test
dds <-DESeq(dds)

