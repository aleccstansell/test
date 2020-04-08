##Set base directories 
work_dir <- "~/doorknobdave/alecS/raw/Healthy_Cells/comparison"
setwd(work_dir)
meta_dir <- "~/doorknobdave/alecS/raw/Healthy_Cells/meta"

###############################################################
#Sample to covariate construction
###############################################################

#Read in Metadata file for all 30 samples used in this experiment
meta_data <- read.csv(file.path(meta_dir, "healthy_meta.csv"), header = TRUE)
s2c <- meta_data

files <- file.path(work_dir, meta_data$sample_title, "run_info.json")

#Initialise data frame
dataframe <- NULL
library(rjson)
i <- 1
#Select all 30 json files read in and rbind
for(i in 1:9){
  
  JsonData <- fromJSON(file= files[i] )
  dataframe <- rbind(JsonData, dataframe)
  
  
}

#Use sample ids as rownames
rownames(dataframe) <- meta_data$sample

#Write to csv file
outdir <- "~/doorknobdave/alecS/alec_code/kallisto/results"
write.csv(dataframe, file.path(outdir, "kallisto_healthy_cells_json_output.csv"))


