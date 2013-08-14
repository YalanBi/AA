#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 25-07-2013
# first written: 25-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#********************************************* this is the final version for analyzing results of testing skipping exons ^_^ *********************************************#
#**************************************************************** testing algorithm: Wilcox.test / ANOVA! ****************************************************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for skipping exons
seMatrix <- NULL
for(chr in 1:5){
  if(file.exists(paste0("Data/AS/SE_chr", chr, "_wt_less_p3.txt"))){
    seMatrix <- rbind(seMatrix, read.table(paste0("Data/AS/SE_chr", chr, "_wt_less_p3.txt"), row.names=NULL))
  } else cat("chr", chr, "NO test!\n")
}
#Bonferroni correction
nrow(seMatrix)# = 5688 exons were tested
length(unique(seMatrix[,1]))# = 4004 genes were tested
testableGenes=unique(seMatrix[,1])
save(testableGenes, file="Data/AS/SE_testableGenes.Rdata")


#calculate the numbers of sig exons and genes
seThre=46.02
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to seThre in each env and across envs from chr1-chr5
matrixGENE <- NULL #a matrix for numbers of genes having at least one exon that -log10(P) are higher than or equal to seThre in each env and across envs from chr1-chr5
for(chr in 1:5){
  if(file.exists(paste0("Data/AS/SE_chr", chr, "_wt_less_p3.txt"))){
    sechr <- read.table(paste0("Data/AS/SE_chr", chr, "_wt_less_p3.txt"), row.names=NULL)
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
    if(length(gn_mixEnv) > 0) save(seGeneList, file=paste0("Data/AS/SE_chr", chr, "_wt_p3.Rdata"))
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
