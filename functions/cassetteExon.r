#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 02-05-2013
# first written: 02-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
#load environment file
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

#cassette exon test
testCassette <- function(expdata = rawexp[exonID, ind_env + 16], ind, use){
  exonprobe <- apply(expdata[ind, ], 2, use)
  #cat(exonprobe, "\n")
  otherprobe <- apply(expdata[!ind, ], 2, use)
  #cat(otherprobe, "\n")
  return(-log10(t.test(exonprobe, otherprobe, alternative="less")$p.value))
}


#*************************************************************** load part ***************************************************************#
#load exp genes
load(file="Data/fullModeMapping/expGenes.Rdata")


#mean or median, change!!!
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  rownameList <- NULL
  
  for(filename in genenames){
    cat(filename, "...\n")
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename), row.names=1, header=T)
    #cat("rawexp loading succeed!\n")
    
    #get probes' and exons' ID
    probes_dir <- probesDir(rawexp)
    #cat("probes of right direction:", probes_dir, "\n")
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
    #cat("exons of right direction:", exonID, "\n")
    
    #uniqueExon <- all tu names of exon probes
    uniqueExon <- unique(rawexp[exonID, "tu"])
    #cat("tu names:", as.character(uniqueExon), "\n")
    
    if(length(uniqueExon) > 1){
      rownameList <- c(rownameList, paste0(gsub(".txt", "", filename), "_", uniqueExon))
      for(exon in uniqueExon){
        #ind <- judge which probe in exonID is of current exon name (T/F)
        ind <- rawexp[exonID, "tu"] == exon
        #cat(as.character(exon), "has probes", exonID[ind], "\n")
        
        res <- NULL
        for(env in 1:4){
          ind_env <- which(as.numeric(menvironment) == env)
          
          #use mean/median to test cassette
          res <- c(res, testCassette(expdata = rawexp[exonID, ind_env + 16], ind, use = mean))
        }
         resmatrix <- rbind(resmatrix, res)
      }
    }
  }
  rownames(resmatrix) <- rownameList
  colnames(resmatrix) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
  write.table(resmatrix, file=paste0("Data/cassetteExon/cassetteExon_chr", chr, "_mean.txt"), sep="\t") #change!!! 
  
  et <- proc.time()[3]
  cat("and finished in", et-st, "s\n")
}
