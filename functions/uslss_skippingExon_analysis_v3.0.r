#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 12-07-2013
# first written: 28-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#********************************************* this is the final version for analyzing results of testing skipping exons ^_^ *********************************************#
#**************************************************************** testing algorithm: Wilcox.test / ANOVA! ****************************************************************#

setwd("D:/Arabidopsis Arrays")

whichFile="_Ex+In"
#calculate the threshold for skipping exons
seMatrix <- NULL
for(chr in 1:5){
  if(file.exists(paste0("Data/AS/SE_chr", chr, whichFile, ".txt"))){
    seMatrix <- rbind(seMatrix, read.table(paste0("Data/AS/SE_chr", chr, whichFile, ".txt"), row.names=NULL))
  } else cat("chr", chr, "NO test!\n")
}
#Bonferroni correction
nrow(seMatrix)# = 4188 exons were tested
-log10(0.05/nrow(seMatrix)/8)# = 5.83; 4188 exons were tested * 4 Env * 2 tests; => seThre=5.83
length(unique(seMatrix[,1]))# = 2761 genes were tested

#calculate the numbers of sig exons and genes
seThre=round(-log10(0.05/nrow(seMatrix)/8), digits=2)# =5.83
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to seThre in each env and across envs from chr1-chr5
matrixGENE <- NULL #a matrix for numbers of genes having at least one exon that -log10(P) are higher than or equal to seThre in each env and across envs from chr1-chr5
for(chr in 1:5){
  if(file.exists(paste0("Data/AS/SE_chr", chr, whichFile, ".txt"))){
    sechr <- read.table(paste0("Data/AS/SE_chr", chr, whichFile, ".txt"), row.names=NULL)
    seGeneList <- list()
    resmatrix <- NULL
    for(e in 1:nrow(sechr)){
      res <- NULL
      for(env in 1:4){
        if(sechr[e, 2*env+1] < seThre && sechr[e, 2*env+2] >= seThre) res <- c(res, 1)
        else res <- c(res, 0)
      }
      resmatrix <- rbind(resmatrix, res)
    }
    
    nTU <- NULL #number of exons that -log10(P) are higher than or equal to seThre in each env
    nGENE <- NULL #number of genes having one or more exons that -log10(P) are higher than or equal to seThre in each env
    for(env in 1:4){
      nTU <- c(nTU, sum(resmatrix[ ,env]))
      seGeneList[[paste0("env", env)]] <- unique(sechr[resmatrix[ ,env] > 0, 1]) #genes that have one or more exons that are spliced out in each env
      nGENE <- c(nGENE, length(seGeneList[[paste0("env", env)]]))
    }
    
    #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to seThre in ANY env
    cnt_mixEnv <- 0
    #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to seThre in ANY env
    gn_mixEnv <- NULL
    for(e in 1:nrow(sechr)){
      if(any(resmatrix[e, ] > 0)){
        cnt_mixEnv <- cnt_mixEnv+1
        gn_mixEnv <- c(gn_mixEnv, sechr[e, 1])
      }
    }
    seGeneList$mixEnv <- unique(gn_mixEnv)
    
    matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
    matrixGENE <- rbind(matrixGENE, c(nGENE, length(seGeneList$mixEnv)))
    if(length(gn_mixEnv) > 0) save(seGeneList, file=paste0("Data/AS/SE_chr", chr, whichFile, ".Rdata"))
  } else{
    cat("chr", chr, "NO test!\n")
    matrixTU <- rbind(matrixTU, c(0,0,0,0,0))
    matrixGENE <- rbind(matrixGENE, c(0,0,0,0,0))
  }
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU

rownames(matrixGENE) <- paste0("chr", 1:5)
colnames(matrixGENE) <- c(paste0("Env", 1:4), "mixEnv")
matrixGENE
