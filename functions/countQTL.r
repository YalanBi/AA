#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 15-04-2013
# first written: 15-04-2013
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


#*************************************************************** load part ***************************************************************#
load(file="Data/fullModeMapping/expGenes.Rdata")


#************************************************************* test and count *************************************************************#
#count main eQTL(qtl >= 4 && int < 5)
countMainQTL <- function(chr = 1, expGeneList, threshold_qtl = 4, threshold_int = 5, cutoffratio = 0.6){
  location <- paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  genenames <- expGeneList[[chr]]
  
  nGene <- 0
  #filename="AT1G01010.txt"
  for(filename in genenames){
    #cat(filename, "starts...\n")
    #load rawexp file, qtl file and int file
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename), row.names=1, header=T)
    #cat("rawexp loading succeed!\n")
    qtl <- read.table(paste0(location, gsub(".txt", "_FM_QTL.txt", filename)), row.names=1, header=T)
    #cat("qtl loading succeed!\n")
    int <- read.table(paste0(location, gsub(".txt", "_FM_Int.txt", filename)), row.names=1, header=T)
    #cat("int loading succeed!\n")
    
    #get probes' and exons' ID
    probes_dir <- probesDir(rawexp)
    #cat("probes of right direction:", probes_dir, "\n")
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
    #cat("exons of right direction:", exonID, "\n")
    
    #uniqueExon <- all tu names of exon probes
    uniqueExon <- unique(rawexp[exonID, "tu"])
    #cat("tu names:", as.character(uniqueExon), "\n")
    
    #to continue testing, or skip and report this gene
    continue <- TRUE
    
    cnt_tu <- 1
    while(cnt_tu %in% 1:length(uniqueExon) && continue){
      #tuID <- ID of exon probes of current tu name
      tuID <- exonID[rawexp[exonID, "tu"] == uniqueExon[cnt_tu]]
      #cat(as.character(uniqueExon[cnt_tu]), "has probes", tuID, "\n")
      
      #if there are at least 3 probes in this tu, then we test
      #if(length(tuID) >= 3){
        #cat(as.character(uniqueExon[cnt_tu]), "has >= 3 probes!\n")
        
        m <- 1
        while(m %in% 1:716 && continue){
          
          #threshold_qtl=4, threshold_int=5, cutoffratio=60%
          #if there are at least 60% probes of this tu are qtl>=4 and int<5, then we stop finding within current gene, and nGene+1
          if(length(which(qtl[tuID, m] >= threshold_qtl && int[tuID, m] < threshold_int))/length(tuID) >= cutoffratio){
            cat(filename, as.character(uniqueExon[cnt_tu]), "marker", m, "has main eQTL!!! and quit from", filename, "!\n")
            continue <- FALSE
            nGene <- nGene + 1
          } else{
            m <- m+1
          }
        }
      #}
      cnt_tu <- cnt_tu + 1
    }
  }
  et <- proc.time()[3]
  cat("chr", chr, "has", nGene, "main eQTL, and finished in", et-st, "s!\n")
  
  return(nGene)
}

mainList <- NULL
for(chr in 1:5){
  mainList <- c(mainList, countMainQTL(chr, expGeneList, threshold_qtl = 4, threshold_int = 5, cutoffratio = 0.6))
}
mainList

#count interaction eQTL(int >= 5)
countInt <- function(chr = 1, expGeneList, threshold_int = 5, cutoffratio = 0.6){
  location <- paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  genenames <- expGeneList[[chr]]
  
  nGene <- 0
  #filename="AT1G01010.txt"
  for(filename in genenames){
    #cat(filename, "starts...\n")
    #load rawexp file and int file
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename), row.names=1, header=T)
    #cat("rawexp loading succeed!\n")
    int <- read.table(paste0(location, gsub(".txt", "_FM_Int.txt", filename)), row.names=1, header=T)
    #cat("int loading succeed!\n")
    
    #get probes' and exons' ID
    probes_dir <- probesDir(rawexp)
    #cat("probes of right direction:", probes_dir, "\n")
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
    #cat("exons of right direction:", exonID, "\n")
    
    #uniqueExon <- all tu names of exon probes
    uniqueExon <- unique(rawexp[exonID, "tu"])
    #cat("tu names:", as.character(uniqueExon), "\n")
    
    #to continue testing, or skip and report this gene
    continue <- TRUE
    
    cnt_tu <- 1
    while(cnt_tu %in% 1:length(uniqueExon) && continue){
      #tuID <- ID of exon probes of current tu name
      tuID <- exonID[rawexp[exonID, "tu"] == uniqueExon[cnt_tu]]
      #cat(as.character(uniqueExon[cnt_tu]), "has probes", tuID, "\n")
      
      #if there are at least 3 probes in this tu, then we test
      #if(length(tuID) >= 3){
        #cat(as.character(uniqueExon[cnt_tu]), "has >= 3 probes!\n")
        
        m <- 1
        while(m %in% 1:716 && continue){
          
          #threshold_int=5, cutoffratio=60%
          #if there are at least 60% probes of this tu are int<5, then we stop finding within current gene, and nGene+1
          if(length(which(int[tuID, m] >= threshold_int))/length(tuID) >= cutoffratio){
            cat(filename, as.character(uniqueExon[cnt_tu]), "marker", m, "has Int!!! and quit from", filename, "!\n")
            continue <- FALSE
            nGene <- nGene + 1
          } else{
            m <- m+1
          }
        }
      #}
      cnt_tu <- cnt_tu + 1
    }
  }
  et <- proc.time()[3]
  cat("chr", chr, "has", nGene, "Int, and finished in", et-st, "s!\n")
  
  return(nGene)
}

intList <- NULL
for(chr in 1:5){
  intList <- c(intList, countInt(chr, expGeneList, threshold_int = 5, cutoffratio = 0.6))
}
