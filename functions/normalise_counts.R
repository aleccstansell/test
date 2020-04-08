#install packages
#install.packages(c("impute","dynamicTreeCut","qvalue","flashClust","Hmisc","WGCNA"))


library(DESeq2)



#Load the DESeq2 package
#suppressPackageStartupMessages(library("DESeq2"))

inputdir <- "~/alec/data/intermediary_data/Kallisto_all"

#Read in count data, with non-unique gene identifiers
countData <- read.csv(file.path(inputdir, "meta_for_coexpression_30.csv"), header = TRUE)
#countData <- read.csv(file.path(inputdir, "meta_for_coexpression_30.csv"), header = TRUE)

#metadir <- "~/alec/data/meta"
#countData <- read.table(file.path(metadir, "metadata_final_30.csv"), header = TRUE, sep = "," row.names = NULL)

#Subset countData data frame, taking only sample counts
geneCounts <- countData[ , 1:30]

#Convert geneCounts data frame to matrix
geneCounts <- as.matrix(geneCounts)

#Convert counts in geneCounts matrix to integers by column
countMatrix <- (apply(geneCounts, MARGIN = 2, FUN = as.integer))

#Subset countData data frame, taking non-unique gene identifiers
genes <- rownames(countData)

#Add unique gene identifiers to countMatrix, as row names
rownames(countMatrix) <- genes

#Make non-unique gene identifiers unique
#genes <- make.unique(genes)

#Read in metadata
metadir <- "~/alec/data/meta"
allMetadata <- read.table(file.path(metadir, "meta_for_coexpression_30.csv"), header = TRUE, sep = ",", col.names = c("sample", "cell", "individual"), row.names = NULL, colClasses = c("character", "factor", "factor"))

#BiocManager::install("DESeq2")
library(DESeq2)
#Create dds in order to run normalisation
dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = allMetadata, design = ~ 1)

#Estimate size factors
dds <- estimateSizeFactors( dds )
sizeFactors(dds)


# dds <- varianceStabilizingTransformation(dds)
# getVarianceStabilizedData(dds)


#Calculate colsums to compare between samples - large differences indicate differences in depth
colSums(counts(dds))

#Plot column sums according to size factor
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))

#The argument normalized equals true, divides each column by its size factor.
logcounts <- log2( counts(dds, normalized=TRUE) + 1 )


## R log transformation

rld <- rlogTransformation(dds)


#Write normalised counts to file
write.table(logcounts, file = "normalised_counts.csv", quote = FALSE, sep = ",") 



logcounts 

plot(logcounts)
pc <- prcomp( t( logcounts ) )

#Remove rows with less than 2 counts per row
dds <- BTdds[rowSums(counts(BTdds)) > 10, ]


dds <- estimateSizeFactors( dds )
sizeFactors(dds)

#Perform DE gene analysis
BTdds <- DESeq(BTdds)

#Get results
BTres <- results(BTdds)




#################




#load packages
library(impute)
library(dynamicTreeCut)
library(qvalue)
library(Hmisc)
library(WGCNA)

#Set wd 
work_dir <- "~/doorknobdave/alecS/WGCNtest/data/metaAnalysisFiles"
setwd(work_dir)


load("metaAnalysisData.RData") 


dataEx
B

B[1] == 2 & B[1,2] == 4


for (i in 1:length(B))
{
  if (B[i]==2)
    print("B is zero")
  else 
    print("B is not zero")
  
}


for (set in 1:exprSize$nSets)
{
  if (sum(!gsg$goodSamples[[set]]))
    printFlush(paste("In set", setLabels[set], "removing samples",
                     paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
  # Remove the offending genes and samples
  multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
}