##  Count number of cell types ##########


#Read metadata in with a column of cell types you want to count
meta_data <- read.csv(file.path(meta_dir, "metadataallqc68.csv"), header = TRUE, stringsAsFactors = F)

#Use the strings as factors = TRUE to get the levels of cell types to see the different cells in dataset
#meta_data <- read.csv(file.path(meta_dir, "metadataallqc68.csv"), header = TRUE, stringsAsFactors = TRUE)
#meta_data$cell

# These are the types in my dataset CM CMPD EM FH LN N NGC TM
  
CM = 19
CMPD = 5
EM = 16
FH = 1
LN = 1
N = 18
NGC = 2
TM = 14

library(stringr)

# Count the number of 'a's in each element of string
numberofcells <- str_count(meta_data$cell, "TM")
sum(numberofcells)

q.data