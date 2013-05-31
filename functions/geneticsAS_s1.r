#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 31-05-2013
# first written: 31-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]


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


ce_threshold = 5.86; threshold_qtl = 8.0; threshold_int = 11.6; cutoffratio = 0.6
for(chr in 1:1){
  location <- paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  cematrix <- read.table(paste0("Data/cassetteExon/cassetteExon_chr", chr, "_allind.txt"), row.names=1, header=T)
  
  #for all genes which have cassette exons in at least one env
  genenames <- unique(unlist(lapply(strsplit(rownames(which(cematrix >= ce_threshold, arr.ind=T)), "_"), "[[", 1)))
  
  genesGAS <- NULL
  genesIAS <- NULL
  
  #filename(AT1G01010)
  for(filename in genenames){
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
    #cat("rawexp loading succeed!\n")
    qtl <- read.table(paste0(location, filename, "_FM_QTL.txt"), row.names=1, header=T)
    #cat("qtl loading succeed!\n")
    int <- read.table(paste0(location, filename, "_FM_Int.txt"), row.names=1, header=T)
    #cat("int loading succeed!\n")
    
    #get probes' and exons' ID
    probes_dir <- probesDir(rawexp)
    #cat("probes of right direction:", probes_dir, "\n")
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
    #cat("exons of right direction:", exonID, "\n")
    
    #uniqueExon <- all tu names of exon probes
    uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
    #cat("tu names:", as.character(uniqueExon), "\n")
    
    for(ce in unlist(lapply(strsplit(grep(filename, rownames(cematrix), value=TRUE), "_"), "[[", 2))){
      tuID <- exonID[rawexp[exonID, "tu"] == ce]
      
      consSigQTLMarkers <- NULL
      consSigIntMarkers <- NULL
      for(m in 1:716){
        if(length(which(qtl[tuID, m] >= threshold_qtl && int[tuID, m] < threshold_int))/length(tuID) >= cutoffratio){
          consSigMarkers <- c(consSigMarkers, m)
        }
        if(length(which(int[tuID, m] >= threshold_int))/length(tuID) >= cutoffratio){
          consSigIntMarkers <- c(consSigIntMarkers, m)
        }
      }
      
      if(length(consSigQTLMarkers) > 0){
        genesGAS <- c(genesGAS, paste0(filename, "_", ce))
      }
      if(length(consSigIntMarkers) > 0){
        genesIAS <- c(genesIAS, paste0(filename, "_", ce))
      }
    }
  }
}


