'Venn Diagram For DESEq2 Results'


#Load Venn Script
source('~/doorknobdave/alecS/alec_code/functions/venn_diagram.R')

#Define Files
file_out_deseq <- "~/doorknobdave/alecS/alec_code/degs_output/HIV_control/deseq2_degs"
load(file.path(file_out_deseq,"final8significantgenes.RData"))


rownames(deseqfinal8signficant$EMvsCM)
setwd('~/doorknobdave/alecS/alec_code/degs_output/HIV_control/venn_diagrams')

vennit(rownames(deseqfinal8signficant$EMvsCM), rownames(deseqfinal8signficant$NvsCM), rownames(deseqfinal8signficant$NvsEM), "EM vs CM", "N vs CM", "N vs EM", "deseq2_comparison_venn.png")

venn

deseqfinal8signficant$EMvsCM