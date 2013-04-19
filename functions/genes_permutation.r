#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 15-04-2013
# first written: 15-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]


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


#*************************************************************** load part ***************************************************************#
load(file="Data/fullModeMapping/expGenes.Rdata")


#************************************************************ make matrix part ************************************************************#
for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  genenames <- expGeneList[[chr]]
  
  newprobematrix <- NULL
  colnameList <- NULL
  
  #filename="AT1G01010.txt"
  for(filename in genenames){
    #cat(filename, "\n")
    rawexp <- read.table(paste0(location, filename), row.names=1, header=T)
    
    if(length(colnameList) == 0){
      colnameList <- c("gene", "tu", "bp", colnames(rawexp)[17:164])
      cat("get colnames!\n")
    }
    
    probes_dir <- probesDir(rawexp)
    exonID <- probes_dir[which(grepl("tu", rawexp[probes_dir, "tu"]))]
    
    newprobe <- filename
    
    for(tu in unique(rawexp[exonID,"tu"]){
      probes4tu <- exonID[which(rawexp[exonID,"tu"] == tu)]
      
      newprobe <- c(newprobe, as.character(tu), round(mean(rawexp[probes4tu, "bp"])))
      
      for(i in 17:164){
        newprobe <- c(newprobe, mean(rawexp[probes4tu, i]))
      }
      newprobematrix <- rbind(newprobematrix, newprobe)
    }
    

  }
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s!\n")
  
  colnames(newprobematrix) <- colnameList
  write.table(newprobematrix, file="Data/fullModeMapping/expGenes_chr", chr, "_1tu1probe.txt", row.names=1, header = TRUE)
}
#to load this file
read.table("Data/fullModeMapping/expGenes_1gene1probe.txt", row.names=1, header=T)
