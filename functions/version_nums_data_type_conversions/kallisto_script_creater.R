

#Read in data from EBI
files_healthy_cells <- read.csv(file.path("~/doorknobdave/alecS/raw/latencystudy","PRJNA369086.txt" ), sep = "\t", header = TRUE)

#write.csv(files_healthy_cells, file.path("~/doorknobdave/alecS/raw/Healthy_Cells","files_report.csv" ))

#Add fastq.gz to files
file1 <- paste(files_healthy_cells$run_accession, "_1.fastq.gz", sep = "")
file2 <- paste(files_healthy_cells$run_accession, "_2.fastq.gz", sep = "")

#Create sample ids that are informative
title <- str_remove_all(files_healthy_cells$sample_title, ' ')



#Make kallisto script
kallisto_script <- paste("kallisto",	"quant",	"-i",	"index.idx",	"-o", title, "-b",	"100", file1, file2, sep = "  ")

#Write to file
setwd("~/doorknobdave/alecS/raw/latencystudy")
write.csv(title, "sample_names.csv")

kallisto  quant  -i  index.idx  -o  Na簿ve-Replicate1  -b  100  SRR4011048_1.fastq.gz  SRR4011048_2.fastq.gz
kallisto  quant  -i  index.idx  -o  CentralMemory-Replicate1  -b  100  SRR4011049_1.fastq.gz  SRR4011049_2.fastq.gz
kallisto  quant  -i  index.idx  -o  EffectorMemory-Replicate1  -b  100  SRR4011050_1.fastq.gz  SRR4011050_2.fastq.gz
kallisto  quant  -i  index.idx  -o  Na簿ve-Replicate2  -b  100  SRR4011051_1.fastq.gz  SRR4011051_2.fastq.gz
kallisto  quant  -i  index.idx  -o  CentralMemory-Replicate2  -b  100  SRR4011052_1.fastq.gz  SRR4011052_2.fastq.gz
kallisto  quant  -i  index.idx  -o  EffectorMemory-Replicate2  -b  100  SRR4011053_1.fastq.gz  SRR4011053_2.fastq.gz
kallisto  quant  -i  index.idx  -o  Na簿ve-Replicate3  -b  100  SRR4011054_1.fastq.gz  SRR4011054_2.fastq.gz
kallisto  quant  -i  index.idx  -o  CentralMemory-Replicate3  -b  100  SRR4011055_1.fastq.gz  SRR4011055_2.fastq.gz
kallisto  quant  -i  index.idx  -o  EffectorMemory-Replicate3  -b  100  SRR4011056_1.fastq.gz  SRR4011056_2.fastq.gz








"kallisto	quant	-i	index.idx	-o	LIR_01-LN_CMPD1lo57lo	-b	100"	







files_healthy_cells$fastq_ftp
strsplit(files_healthy_cells$fastq_ftp,"/")


names <- strsplit(as.character(files_healthy_cells$fastq_ftp), split = "/")

names <- as.data.frame(names)
names[,8]

names[,13]



fruits <- c("one apple", "two pears", "three bananas")
str_remove(fruits, "[aeiou]")
#> [1] "ne apple"     "tw pears"     "thre bananas"
str_remove_all(fruits, "[aeiou]")
#> [1] "n ppl"    "tw prs"   "thr bnns"