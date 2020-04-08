#' Remove version number from ensembl gene/transcript ids
#'
#' @param ids vector of ensembl ids
#' @return vector of ensembl ids without the version number
#' @import stringr
#' @export
#' @examples
#' removeVersion("ENSMUSG00000001017.15")
#' @author Beth Signal
#' 
#' 
removeVersion <- function(ids){
  return(unlist(lapply(stringr::str_split(ids, "[.]"), "[[",1)))
}

library(praise)


##Set base directories 

work_dir <- "~/doorknobdave/alecS/intermediary_data/kallisto_out"
meta_dir <-  "~/doorknobdave/alecS/meta"

#Read in sample meta data for file names

meta <- read.csv(file.path(meta_dir, "metadataallqc68_no_spaces.csv"), header = TRUE)

#Create object file names containing file paths for all the files 
file.names <- file.path(work_dir, meta$sample, "abundance.tsv")

file.location <-  file.path(work_dir, meta$sample)

#ids <- read.table("~/doorknobdave/alecS/intermediary_data/kallisto_out/LIR_01-PBMC_CM/abundance.tsv", sep = "\t", header = TRUE)

# Loop through all files creating a new tsv with no version names

file.names[10]

for(i in 1:length(file.names)){
  ids <- read.table(file.names[i], sep = "\t", header = TRUE)
  ids$target_id <- removeVersion(ids$target_id)
  setwd(file.location[i])
  write.table(ids, "abundance_no_versions.tsv", quote = FALSE, sep = "\t")
  print(paste("In progress....", i/30*100 , "%", "of the way", sep = " ", collapse = NULL))
  print(praise(template = "fricken ${adjective}!"))
  if (i == 30) {
    print(praise("${Exclamation}! You have ${adjective} version-number free transcript ids!"))
  } 
  
}