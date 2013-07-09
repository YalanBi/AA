#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 01-07-2013
# first written: 28-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#****************************** this is the final version for analyzing results of testing GENETIC/INTERACTION regulated skipping exons ^_^ ******************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for INTERACTION regulated skipping exons
iseMatrix <- NULL
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/SE_chr", chr, "_IAS_wt_less.txt"))){
    iseMatrix <- rbind(iseMatrix, read.table(paste0("Data/geneticsAS/SE_chr", chr, "_IAS_wt_less.txt"), row.names=NULL))
  } else cat("chr", chr, "NO test!\n")
}
#Bonferroni correction
-log10(0.05/nrow(iseMatrix)/4)# = 3.35; 28 exons were tested * 4 Env; => iseThre=3.35
rm(iseMatrix)

#calculate the numbers of sig exons and genes
iseThre=3.35
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to iseThre in each env and across envs from chr1-chr5
matrixGENE <- NULL #a matrix for numbers of genes having at least one exon that -log10(P) are higher than or equal to iseThre in each env and across envs from chr1-chr5
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/SE_chr", chr, "_IAS_wt_less.txt"))){
    isechr <- read.table(paste0("Data/geneticsAS/SE_chr", chr, "_IAS_wt_less.txt"), row.names=NULL)
    iseGeneList <- list()
    resmatrix <- NULL
    
    for(e in 1:nrow(isechr)){
      res <- NULL
      for(env in 1:4){
        if((isechr[e,env*2+2] >= iseThre && isechr[e,env*2+3] < iseThre) || (isechr[e,env*2+2] < iseThre && isechr[e,env*2+3] >= iseThre)){
          res <- c(res, 1)
        } else res <- c(res, 0)
      }
      resmatrix <- rbind(resmatrix, res)
    }
    
    nTU <- NULL #number of exons that -log10(P) are higher than or equal to iseThre in one gt and lower in the other gt, each env sep
    nGENE <- NULL #number of genes having one or more exons that -log10(P) are higher than or equal to iseThre in each env
    for(env in 1:4){
      nTU <- c(nTU, sum(resmatrix[ ,env]))
      iseGeneList[[paste0("env", env)]] <- unique(isechr[resmatrix[ ,env] > 0, 1]) #genes that have one or more exons that are spliced out in each env
      nGENE <- c(nGENE, length(iseGeneList[[paste0("env", env)]]))
    }
    
    #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to iseThre in ANY env
    cnt_mixEnv <- 0
    #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to iseThre in ANY env
    gn_mixEnv <- NULL
    for(e in 1:nrow(isechr)){
      if(sum(resmatrix[e, ]) > 0){
        cnt_mixEnv <- cnt_mixEnv + 1
        gn_mixEnv <- c(gn_mixEnv, isechr[e,1])
      }
    }
    iseGeneList$mixEnv <- unique(gn_mixEnv)
    
    matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
    matrixGENE <- rbind(matrixGENE, c(nGENE, length(iseGeneList$mixEnv)))
    save(iseGeneList, file = paste0("Data/geneticsAS/SE_chr", chr, "_IAS_gnList_wt.Rdata"))
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
