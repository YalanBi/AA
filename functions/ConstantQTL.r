#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 12-04-2013
# first written: 12-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)


#direction selection
probesDir <- function(exp_data = rawexp){
  if(unique(exp_data[,"strand"]) == "sense"){
    direction_id <- which(exp_data[, "direction"] == "reverse")
  }
  if(unique(exp_data[,"strand"]) == "complement"){
    direction_id <- which(exp_data[, "direction"] == "forward")
  }
  return(direction_id)
}





findExpGene <- function(chr, threshold = 5){
  location <- paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  expGene <- NULL
  genenames <- dir(location)[which(grepl("_QTL.txt", dir(location))]
  
  #filename "AT1G01010_FM_QTL.txt"
  for(filename in genenames){
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", gsub("_FM_QTL", "", filename)), row.names=1, header=T)
    
    probes_dir <- probesDir(rawexp)
    #cat(filename, "\nprobeDir:", probes_dir, "\n")
    #exonID <- exons of right direction
    exonID <- probes_dir[which(grepl("tu", rawexp[probes_dir, "tu"]))]
    #cat("exons:", exonID, "\n")
    #intronID <- introns of right direction
    intronID <- probes_dir[which(grepl("intron", rawexp[probes_dir, "tu"]))]
    #cat("introns:", intronID, "\n")
    
    if(length(exonID) > 0){
      cat(filename, "median is", median(unlist(rawexp[exonID, 17:164])))
      if(median(unlist(rawexp[exonID, 17:164])) >= threshold){
        expGene <- c(expGene, filename)
        cat(" TRUE!")
      }
      cat("\n")
    }
  }
  et <- proc.time()[3]
  cat("finding expressed genes on chr", chr, "finished in", et-st, "s\n")
  
  return(expGene)
}
