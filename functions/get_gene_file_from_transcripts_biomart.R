
#Function to return list of mapping genes to transcripts
#May be used in tximport for sumarising to gene level or in sleuth for determining gene level abundance

#Load Biomart Package
library(biomaRt)

#No inputs required returns biomart mappings for all transcripts within hsapiens
gett2g <- function(){
  
  #Download mart object for use
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "hsapiens_gene_ensembl")
  
  #Select relevant information needed 
  t2g <- biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "transcript_version",
                   "ensembl_gene_id", "external_gene_name", "description",
                   "transcript_biotype"),
    mart = mart)
  
  #Rename columns
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  
  #Pull out relevant ones of interest
  t2g <- dplyr::select(t2g, c('target_id', 'ens_gene', 'ext_gene'))
  
  
  return(t2g)
  
  
  
}


