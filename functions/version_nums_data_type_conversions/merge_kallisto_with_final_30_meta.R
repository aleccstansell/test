dir <- "~/doorknobdave/alecS/alec_code/kallisto"

setwd(dir)

meta <- read.csv("metadata_final_30.csv")
meta$

kal <- NULL
kal <- read.csv("kallistoscriptall.csv", header = F)
head(kal)

kal


getwd()

merger <- merge(meta, kal, by.x="sample", by.y ="V6")

head(merger)
write.csv(merger, "kallisto_script_30.csv")
