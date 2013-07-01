#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 01-07-2013
# first written: 19-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#************************************************* this is the final version for summarizing probes into exon level! ^_^ *************************************************#
#************************************** one probe for one exon, and the order of RILs is as same as their initial order (no change) **************************************#
#********************************************* this is used for PERMUTAION, CIS-TRANS PLOT and MAIN/CONSISTENT eQTL and Int! *********************************************#

setwd("D:/Arabidopsis Arrays")
menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]
#load exp genes
load(file="Data/ExpGenes/expGenes_final.Rdata")

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

#summarize probes in one exon into one probe, each exon has one probe only, chr1-5 separately
#make new expression files among expressed genes, 1 tu 1 probe, at exon level
for(chr in 1:5){
  location <- paste0("Data/Raw/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  genenames <- expGeneList[[chr]]
  
  newprobematrix <- NULL
  rownameList <- NULL
  colnameList <- NULL
  
  #filename="AT1G01010"
  for(filename in genenames){
    rawexp <- read.table(paste0(location, filename, ".txt"), row.names=1, header=TRUE)
    
    if(length(colnameList) == 0){
      colnameList <- c("bp", colnames(rawexp)[17:164])
      cat("get colnames!\n")
    }
    
    probes_dir <- probesDir(rawexp)
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
    uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
    
    #cat(filename)
    for(tu in uniqueExon){
      probes4tu <- exonID[which(rawexp[exonID, "tu"] == tu)]
      if(length(probes4tu) > 0){
        rownameList <- c(rownameList, paste0(filename, "_", tu))
        newprobematrix <- rbind(newprobematrix, c(round(mean(rawexp[probes4tu, "bp"])), colMeans(rawexp[probes4tu, 17:164])))
        #cat(" ", tu)
      } #else cat(paste0(" NO", tu))
    }
    #cat("\n")
  }
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s!\n\n")
  
  rownames(newprobematrix) <- rownameList
  colnames(newprobematrix) <- colnameList
  write.table(newprobematrix, file=paste0("Data/summarizedGene/expGenes_chr", chr, "_1tu1probe.txt"), sep="\t")
}

#to load this file
read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_1tu1probe.txt"), row.names=1, header=TRUE)
