###################################################################################################################################
#START




##Set base directories 

work_dir <- "~/doorknobdave/alecS/raw/Healthy_Cells/abundance"
meta_dir <-  "~/doorknobdave/alecS/raw/Healthy_Cells/meta"

#Read in sample meta data for file names

s2c <- read.csv(file.path(meta_dir, "healthy_meta.csv"), header = TRUE)

#Create object file names containing file paths for all the files 
file.names <- file.path(work_dir, s2c$sample_title, "abundance.tsv")

all(file.exists(file.names))



############ Might be more effective to use entrez ids the matrix of all counts ...


#### Code to determine whether the genes/transcripts are the same for all the kallisto outputs they should be in the same order as well in order to match the correct genes

#initialise object d
d <- ""
for(i in 1:length(file.names)){
  test <- read.table(file.names[i],header=TRUE, sep="\t", stringsAsFactors=FALSE)
  testb <- read.table(file.names[i+1],header=TRUE, sep="\t", stringsAsFactors=FALSE)
  d <- rbind(d, all(testb$target_id == test$target_id))
  
}

####Code to add all matrices together + to lable each column name as sample name

#Initialise object to use in for loop
mastercol <- NULL

#For loop through all files adding only the estimated counts to a matrix containing ecount data from all the samples
for(i in 1:length(file.names)){
    counts <- read.table(file.names[i],header=TRUE, sep="\t", stringsAsFactors=FALSE) 
    mastercol <- cbind(mastercol, counts$est_counts) # can add counts$target_id to do a visual comparison if all are the same above code checks this though)
    print(praise("${Exclamation}! ${adjective}"))
    
}


#For loop through all files adding only the estimated counts to a matrix containing ecount data from all the samples
# for(i in 1:length(file.names)){
#   counts <- read.table(file.names[i],header=TRUE, sep="\t", stringsAsFactors=FALSE) 
#   mastercol <- cbind(mastercol, counts$est_counts) # can add counts$target_id to do a visual comparison if all are the same above code checks this though)
#   
# }


#Read in data frame of one of the results in order to get the target ids to be used as row names in the final object
names <-read.table(file.names[15],header=TRUE, sep="\t", stringsAsFactors=FALSE)


#set colnames of mastercol to the sample names 
colnames(mastercol) <- s2c$sample

#make matrix
mastercol <- as.matrix(mastercol)

#add rownames for gene symbols

rownames(mastercol) <- names$target_id


##Write to file
setwd(work_dir)
write.table(mastercol, file = "mega_count_matrix.csv", quote = FALSE, sep = ",") 


################### Need code to remove identifiers in R some form of search algorith

####### Basic if statement example to use ########################
a = c(3, 4, 5, 6, .1)
b = 4 

  
  if (a == b) {
    print("a is the same as b")
  } else {
    print("a is not the same as b")
  }

for(i in 1:10){
  if (a[i] == ".1")
}

round(enst344.1)




str_remove(test$target_id, "[.i]" )
#To remove different target ids
i = 1

######DOESNT WORK
#Remove addtional identifiers - for loop this!! Not sure how to use i within .i etc
absc <- c("howzit.1", "sup.2")
absc
absce <- str_remove(absc, ".1")
absce
absce <- str_remove(absce, ".2")
absce



test <- read.table(file.names[1],header=TRUE, sep="\t", stringsAsFactors=FALSE)

write.table(test, file = "test", quote = FALSE, sep = ",")


# target <- str_remove(test$target_id, ".1" )
# head(target)
# target<- str_remove(test$target_id, ".2" )
# target<- str_remove(test$target_id, ".3" )
# target <- str_remove(test$target_id, ".4" )
# target <- str_remove(test$target_id, ".5" )
# target <- str_remove(test$target_id, ".6" )
# target <- str_remove(test$target_id, ".7" )
# target <- str_remove(test$target_id, ".8" )
# target <- str_remove(test$target_id, ".9" )
# target <- str_remove(test$target_id, ".0" )

target
test$target_id

#Remove addtional identifiers - for loop this!! Not sure how to use i within .i etc
test$target_id <- str_remove(test$target_id, "_1]" )
test$target_id <- str_remove(test$target_id, "_2]" )
test$target_id <- str_remove(test$target_id, "_3]" )
test$target_id <- str_remove(test$target_id, "_4]" )
test$target_id <- str_remove(test$target_id, "_5]" )
test$target_id <- str_remove(test$target_id, "_6]" )
test$target_id <- str_remove(test$target_id, "_7]" )
test$target_id <- str_remove(test$target_id, "_8]" )
test$target_id <- str_remove(test$target_id, "_9]" )
test$target_id <- str_remove(test$target_id, "[_0]" )





colnames(mastercol) <-s2c$sample

colnames(mastercol) <- s2c$sample[1:length(colnames(mastercol))]


head(test)




test



head(mastercol)

fix(mastercol)

head(mastercol)

######################################################################################################################################




