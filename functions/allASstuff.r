#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 11-07-2013
# first written: 11-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


setwd("D:/Arabidopsis Arrays")

#all AS stuff by wilcox.test
allASGenes <- list()
for(chr in 1:5){
  gn <- NULL
  load(file=paste0("Data/AS/SE_chr", chr, "_gnList_wt.Rdata"))
  gn <- c(gn, seGeneList$mixEnv)
  rm(seGeneList)
  
  load(file=paste0("Data/AS/5'AS_chr", chr, "_gnList_wt.Rdata"))
  gn <- c(gn, asGeneList$mixEnv)
  rm(asGeneList)
  
  load(file=paste0("Data/AS/3'AS_chr", chr, "_gnList_wt.Rdata"))
  gn <- c(gn, asGeneList$mixEnv)
  rm(asGeneList)
  
  load(file=paste0("Data/AS/RI_chr", chr, "_gnList_wt.Rdata"))
  gn <- c(gn, riGeneList$mixEnv)
  rm(riGeneList)
  
  allASGenes[[paste0("chr", chr)]] <- sort(unique(gn))
}
save(allASGenes, file="Data/AS/allASGenes_wt.Rdata")

#all AS stuff by ANOVA
allASGenes <- list()
for(chr in 1:5){
  gn <- NULL
  load(file=paste0("Data/AS/SE_chr", chr, "_gnList_ANOVA.Rdata"))
  gn <- c(gn, seGeneList$mixEnv)
  rm(seGeneList)
  
  load(file=paste0("Data/AS/5'AS_chr", chr, "_gnList_ANOVA.Rdata"))
  gn <- c(gn, asGeneList$mixEnv)
  rm(asGeneList)
  
  load(file=paste0("Data/AS/3'AS_chr", chr, "_gnList_ANOVA.Rdata"))
  gn <- c(gn, asGeneList$mixEnv)
  rm(asGeneList)
  
  load(file=paste0("Data/AS/RI_chr", chr, "_gnList_ANOVA.Rdata"))
  gn <- c(gn, riGeneList$mixEnv)
  rm(riGeneList)
  
  allASGenes[[paste0("chr", chr)]] <- sort(unique(gn))
}
save(allASGenes, file="Data/AS/allASGenes_ANOVA.Rdata")
