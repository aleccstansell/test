##Transcript abundance to Gene counts (ensembl)

# The following code takes transcript count data outputed from kallisto and summarises the abundance to the gene level for further use in coexpression networks

## Install Packages

#tximport
# try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("tximport")
#rdhf5
# bioccLite("rhdf5")

# Load libraries

library(rhdf5)
library(tximport)

###############WORKING TEST###############

##Set base directories 
work_dir <- "~/doorknobdave/alecS/intermediary_data/kallisto_out"
meta_dir <-  "~/doorknobdave/alecS/meta"


##get transcript to gene id table from ensembl for anotations between genes and transcripts
#http://www.ensembl.org/biomart/martview/af949c4506992d3c4c752b186d6cba44


#Read in tx2 gene file for conversion between transcripts and gene ids
tx2gene <- read.csv(file.path(metadir, "tx2gene94.txt"))
head(tx2gene)

#Directory where sequence files meta data etc are located
#dir <- "~/doorknobdave/alecS/intermediary_data/kallisto_out"

#Read all meta data in table consisting of everything known about samples as an object s2c - sample to covariate matrix 
samples <- read.csv(file.path(meta_dir, "metadata_final_30.csv"), header = TRUE)

#Specify file paths for each of the samples h5 abundance estimate files using the sample title names (the same as the output file names created in kallisto)
files <- file.path(work_dir, samples$sample, "abundance.h5") #abundance.h5


#specify names of each sample
names(files) <- samples$sample

#Check all files exist
all(file.exists(files))

#Create txi object for the conversion
txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)

#txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
tail(txi.kallisto$counts)

#See number of genes transcripts are summarised to 
dim(txi.kallisto$counts)

#Write gene level abundance to file
out_dir <- "~/doorknobdave/alecS/intermediary_data/Kallisto_all"
setwd(out_dir)
write.table(txi.kallisto$counts, "genelevelabundance_30.csv", quote = FALSE, sep = ",")


