#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 03-07-2013
# first written: 03-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#****************************************** this is the final version for summizing genes having GENETIC regulated AS events! ^_^ ******************************************#

setwd("D:/Arabidopsis Arrays")

gasGenes <- list()
for(chr in 1:5){
  #start with skipping exons
  if(file.exists(paste0("Data/geneticsAS/SE_chr", chr, "_GAS_genenameList.Rdata"))){
    load(file = paste0("Data/geneticsAS/SE_chr", chr, "_GAS_genenameList.Rdata"))
    gasGenes$se[[paste0("chr", chr)]] <- gseGeneList[[5]]
    rm(gseGeneList)
    cat("chr", chr, "finished, GENETICS skipping exons!\n")
  } else{
    gasGenes$se[[paste0("chr", chr)]] <- NULL
    cat("chr", chr, "NO genes, GENETICS skipping exons!\n")
  }
  
  #next is 5'AS
  if(file.exists(paste0("Data/geneticsAS/5'AS_chr", chr, "_GAS_genenameList.Rdata"))){
    load(file = paste0("Data/geneticsAS/5'AS_chr", chr, "_GAS_genenameList.Rdata"))
    gasGenes$as5[[paste0("chr", chr)]] <- g5asGeneList[[5]]
    rm(g5asGeneList)
    cat("chr", chr, "finished, GENETICS 5'site AS!\n")
  } else{
    gasGenes$as5[[paste0("chr", chr)]] <- NULL
    cat("chr", chr, "NO genes, GENETICS 5'site AS!\n")
  }
  
  #then is 3'AS
  if(file.exists(paste0("Data/geneticsAS/3'AS_chr", chr, "_GAS_genenameList.Rdata"))){
    load(file = paste0("Data/geneticsAS/3'AS_chr", chr, "_GAS_genenameList.Rdata"))
    gasGenes$as3[[paste0("chr", chr)]] <- g3asGeneList[[5]]
    rm(g3asGeneList)
    cat("chr", chr, "finished, GENETICS 3'site AS!\n")
  } else{
    gasGenes$as3[[paste0("chr", chr)]] <- NULL
    cat("chr", chr, "NO genes, GENETICS 3'site AS!\n")
  }
  
  #end up with RI
  if(file.exists(paste0("Data/geneticsAS/RI_chr", chr, "_GAS_genenameList.Rdata"))){
    load(file = paste0("Data/geneticsAS/RI_chr", chr, "_GAS_genenameList.Rdata"))
    gasGenes$ri[[paste0("chr", chr)]] <- griGeneList[[5]]
    rm(griGeneList)
    cat("chr", chr, "finished, GENETICS retained introns!\n\n")
  } else{
    gasGenes$ri[[paste0("chr", chr)]] <- NULL
    cat("chr", chr, "NO genes, GENETICS retained introns!\n\n")
  }
}
save(gasGenes, file="Data/geneticsAS/AS_byChr_genenameList.Rdata")
