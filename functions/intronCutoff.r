#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 12-04-2013
# first written: 12-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")


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

intronExp <- function(rawexp, intronID, threshould){
  for(p in intronID){
    for(i in 17:164){
      if(rawexp[p, i] >= threshould){
        rawexp[p, i] <- 1
      } else{
        rawexp[p, i] <- 0
      }
    }
  }
  return(rawexp[intronID, 17:164])
}


threshould = 5

nGene <- NULL
nInall <- 0
nInHigh <- 0

for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  cat("\nNow chr", chr, "starts!\n")
  st <- proc.time()[3]
  counting <- 0
  
  #filename "AT1G01010.txt"
  for(filename in dir(location)[which(grepl(".txt", dir(location)) & !grepl("_QTL", dir(location)))]){
    rawexp <- read.table(paste0(location, filename), row.names=1, header=T)
    counting <- counting + 1
    
    probes_dir <- probesDir(rawexp)
    #cat(filename, "\nprobeDir:", probes_dir, "\n")
    #intronID <- introns of right direction
    intronID <- probes_dir[which(grepl("intron", rawexp[probes_dir, "tu"]))]
    #cat("introns:", intronID, "\n")
        
    if(length(intronID) > 0){
      #number of intron individuals higher than threshould
      nInHigh <- nInHigh + length(which(rawexp[intronID, 17:164] >= threshould))
      #number of all intron individuals
      nInall <- nInall + length(intronID)*148
    }
    
    #cat(counting, "\t")
  }
  et <- proc.time()[3]
  cat("chr", chr, "has", counting, "genes; and ends up in", et-st, "s\n")
  
  cat("until chr", chr, "there are", nInHigh, "intron individuals higher than", threshould, "and", nInall, "intron individuals in all\n")
  nGene <- c(nGene, counting)
}



hist(exonmean, breaks=50)
hist(intronmean, breaks=50, col=6, add=TRUE)
