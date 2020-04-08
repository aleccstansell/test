####
#Venn Diagram Code
###

# install.packages('VennDiagram')

# Load library
library(VennDiagram)


# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Dark2")

# Chart

?brewer.pal

Emvs









vennit <- function(set1, set2, set3, nameset1, nameset2, nameset3, name_venn){
  
  
  # Prepare a palette of 3 colors with R colorbrewer:
  library(RColorBrewer)
  myCol <- brewer.pal(3, "Paired")
  
  venn.diagram(
    x = list(set1, set2, set3),
    category.names = c(nameset1, nameset2, nameset3),
    filename = name_venn,
    output=TRUE,
    
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
  )
  
  
}




# # Generate 3 sets of 200 words
# set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# 
# # Generate 3 sets of 200 words
# set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")


# setwd("~/doorknobdave/alecS/alec_code/degs_output")

# vennit(set1, set2,set3,"set1", "set2","set3", "lekkerven3n.png")

# set1 <- 
# 
# set2 <- 
#   
# set3 <- 
# 
# # 
# 
# venn.diagram(
#   x = list(set1, set2, set3),
#   category.names = c("EM_CM_sleuth" , "EM_CM_deseq2" , "EM_CM_limma"),
#   filename = '#14_venn_diagramm.png',
#   output=TRUE,
#   
#   # Output features
#   imagetype="png" ,
#   height = 480 , 
#   width = 480 , 
#   resolution = 300,
#   compression = "lzw",
#   
#   # Circles
#   lwd = 2,
#   lty = 'blank',
#   fill = myCol,
#   
#   # Numbers
#   cex = .6,
#   fontface = "bold",
#   fontfamily = "sans",
#   
#   # Set names
#   cat.cex = 0.6,
#   cat.fontface = "bold",
#   cat.default.pos = "outer",
#   cat.pos = c(-27, 27, 135),
#   cat.dist = c(0.055, 0.055, 0.085),
#   cat.fontfamily = "sans",
#   rotation = 1
# )
# 
