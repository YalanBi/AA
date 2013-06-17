#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 16-06-2013
# first written: 16-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
#load SNP probes' number, snp is a list of 475009 integers
snp <- read.table("snp ids.txt")[,1]


removeSNP <- function(filename = "AT1G01010", location, verbose = FALSE){
  rawexp <- read.table(paste0(location, filename, ".txt"), row.names=1, header=T)

  if(any(rownames(rawexp) %in% snp)){
    if(verbose) cat(filename, "has SNP in it!")

    #save old exp file with name "AT1G01010_SNPin.txt"
    write.table(rawexp, file = paste0(location, filename, "_SNPin.txt"))

    #save new exp file with name "AT1G01010.txt"
    rawexp <- rawexp[-which(rownames(rawexp) %in% snp), ]
    write.table(rawexp, file = paste0(location, filename, ".txt"))

    if(verbose) cat(" and old exp+new exp files saved!\n")
  } else if(verbose) cat(filename, "no SNP\n")
}


#*************************************************************** remove part ***************************************************************#
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts now!\n")
  
  location <- paste0("Data/Raw/chr", chr, "_norm_hf_cor/")
  genenames <- gsub(".txt", "", dir(location))
  
  #filename = "AT1G01010"
  for(filename in genenames){
    removeSNP(filename, location)
  }
  et <- proc.time()[3]
  cat("and chr", chr, "finished in", et-st, "s!\n")
}

