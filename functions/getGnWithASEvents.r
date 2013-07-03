#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 03-07-2013
# first written: 03-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#*************************************************** this is the final version for summizing genes having AS events! ^_^ ***************************************************#
#***************************************************************** AS events: no matter regulated by what! *****************************************************************#

setwd("D:/Arabidopsis Arrays")

asGenes <- list()
for(chr in 1:5){
  #start with skipping exons
  if(file.exists(paste0("Data/AS/SE_chr", chr, "_genenameList.Rdata"))){
    load(file = paste0("Data/AS/SE_chr", chr, "_genenameList.Rdata"))
    asGenes$se[[paste0("chr", chr)]] <- seGeneList[[5]]
    rm(seGeneList)
    cat("chr", chr, "finished, skipping exons!\n")
  } else{
    asGenes$se[[paste0("chr", chr)]] <- NULL
    cat("chr", chr, "NO genes, skipping exons!\n")
  }
  
  #next is 5'AS
  if(file.exists(paste0("Data/AS/5AS_chr", chr, "_genenameList.Rdata"))){
    load(file = paste0("Data/AS/5AS_chr", chr, "_genenameList.Rdata"))
    asGenes$as5[[paste0("chr", chr)]] <- as5GeneList[[5]]
    rm(as5GeneList)
    cat("chr", chr, "finished, 5'site AS!\n")
  } else{
    asGenes$as5[[paste0("chr", chr)]] <- NULL
    cat("chr", chr, "NO genes, 5'site AS!\n")
  }
  
  #then is 3'AS
  if(file.exists(paste0("Data/AS/3AS_chr", chr, "_genenameList.Rdata"))){
    load(file = paste0("Data/AS/3AS_chr", chr, "_genenameList.Rdata"))
    asGenes$as3[[paste0("chr", chr)]] <- as3GeneList[[5]]
    rm(as3GeneList)
    cat("chr", chr, "finished, 3'site AS!\n")
  } else{
    asGenes$as3[[paste0("chr", chr)]] <- NULL
    cat("chr", chr, "NO genes, 3'site AS!\n")
  }

  #end up with RI
  if(file.exists(paste0("Data/AS/RI_chr", chr, "_genenameList.Rdata"))){
    load(file = paste0("Data/AS/RI_chr", chr, "_genenameList.Rdata"))
    asGenes$ri[[paste0("chr", chr)]] <- riGeneList[[5]]
    rm(riGeneList)
    cat("chr", chr, "finished, retained introns!\n\n")
  } else{
    asGenes$ri[[paste0("chr", chr)]] <- NULL
    cat("chr", chr, "NO genes, retained introns!\n")
  }
}
save(asGenes, file="Data/AS/AS_byChr_genenameList.Rdata")
