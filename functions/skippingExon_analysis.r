#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 10-07-2013
# first written: 28-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#********************************************* this is the final version for analyzing results of testing skipping exons ^_^ *********************************************#
#**************************************************************** testing algorithm: Wilcox.test / ANOVA! ****************************************************************#

setwd("D:/Arabidopsis Arrays")

whichTest="wt"# "ANOVA"

#calculate the threshold for skipping exons
seMatrix <- NULL
for(chr in 1:5){
  seMatrix <- rbind(seMatrix, read.table(paste0("Data/AS/SE_chr", chr, "_", whichTest, "_less.txt"), row.names=NULL))
}
#Bonferroni correction
nrow(seMatrix)# = 49813 exons were tested
-log10(0.05/nrow(seMatrix)/4)# = 6.60; 49813 exons were tested * 4 Env; => seThre=6.60
length(unique(seMatrix[,1]))# = 12776 genes were tested
rm(seMatrix)

#calculate the numbers of sig exons and genes
seThre=6.60
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to seThre in each env and across envs from chr1-chr5
matrixGENE <- NULL #a matrix for numbers of genes having at least one exon that -log10(P) are higher than or equal to seThre in each env and across envs from chr1-chr5
for(chr in 1:5){
  if(file.exists(paste0("Data/AS/SE_chr", chr, "_", whichTest, "_less.txt"))){
    sechr <- read.table(paste0("Data/AS/SE_chr", chr, "_", whichTest, "_less.txt"), row.names=NULL)
    seGeneList <- list()
    
    nTU <- NULL #number of exons that -log10(P) are higher than or equal to seThre in each env
    nGENE <- NULL #number of genes having one or more exons that -log10(P) are higher than or equal to seThre in each env
    for(env in 1:4){
      nTU <- c(nTU, length(which(sechr[ ,env+2] >= seThre)))
      seGeneList[[paste0("env", env)]] <- unique(sechr[sechr[ ,env+2] >= seThre, 1]) #genes that have one or more exons that are spliced out in each env
      nGENE <- c(nGENE, length(seGeneList[[paste0("env", env)]]))
    }
    
    #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to seThre in ANY env
    cnt_mixEnv <- 0
    #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to seThre in ANY env
    gn_mixEnv <- NULL
    for(e in 1:nrow(sechr)){
      if(any(sechr[e, 3:6] >= seThre)){
        cnt_mixEnv <- cnt_mixEnv + 1
        gn_mixEnv <- c(gn_mixEnv, sechr[e,1])
      }
    }
    seGeneList$mixEnv <- unique(gn_mixEnv)
    
    matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
    matrixGENE <- rbind(matrixGENE, c(nGENE, length(seGeneList$mixEnv)))
    save(seGeneList, file=paste0("Data/AS/SE_chr", chr, "_gnList_", whichTest, ".Rdata"))
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
