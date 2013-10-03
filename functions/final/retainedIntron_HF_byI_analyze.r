#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 03-09-2013
# first written: 03-09-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#*********************************** this is the final version for analyzing results of testing interaction regulated retained introns ^_^ ***********************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for interaction regulated retained introns
iriMatrix <- NULL
for(chr in 1:5){
  iriMatrix <- rbind(iriMatrix, read.table(paste0("Data/geneticsAS/retainedIntronByI_chr", chr, "_wt_p2.txt"), row.names=NULL))
}
#Bonferroni correction
nrow(iriMatrix)# = 60981 introns were tested 
-log10(0.05/nrow(iriMatrix)/8)# 60981 introns were tested * 8 Env; => iriThre=6.99
length(unique(iriMatrix[,1]))# = 12995 genes were tested

#calculate the numbers of NOT sig introns and genes
iriThre=round(-log10(0.05/82943/8), digits=2)# =7.12
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to iriThre in each env and across envs from chr1-chr5
matrixGENE <- NULL #a matrix for numbers of genes having at least one exon that -log10(P) are higher than or equal to iriThre in each env and across envs from chr1-chr5
for(chr in 1:5){
  irichr <- read.table(paste0("Data/geneticsAS/retainedIntronByI_chr", chr, "_wt_p2.txt"), row.names=NULL)
  iriGeneList <- list()
  resmatrix <- NULL
  
  for(e in 1:nrow(irichr)){
    res <- NULL
    for(env in 1:4){
      if((irichr[e,env*2+2] >= iriThre && irichr[e,env*2+3] < iriThre) || (irichr[e,env*2+2] < iriThre && irichr[e,env*2+3] >= iriThre)){
        res <- c(res, 1)
      }else res <- c(res, 0)
    }
    resmatrix <- rbind(resmatrix, res)
  }
  
  nTU <- NULL #number of exons that -log10(P) are higher than or equal to iriThre in one gt and lower in the other gt, each env sep
  nGENE <- NULL #number of genes having one or more exons that -log10(P) are higher than or equal to iriThre in each env
  for(env in 1:4){
    nTU <- c(nTU, sum(resmatrix[ ,env]))
    iriGeneList[[paste0("env", env)]] <- unique(irichr[resmatrix[ ,env] > 0, 1]) #genes that have one or more exons that are spliced out in each env
    nGENE <- c(nGENE, length(iriGeneList[[paste0("env", env)]]))
  }
  
  #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to iriThre in ANY env
  cnt_mixEnv <- 0
  #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to iriThre in ANY env
  gn_mixEnv <- NULL
  for(e in 1:nrow(irichr)){
    if(sum(resmatrix[e, ]) > 0){
      cnt_mixEnv <- cnt_mixEnv + 1
      gn_mixEnv <- c(gn_mixEnv, irichr[e,1])
    }
  }
  iriGeneList$mixEnv <- unique(gn_mixEnv)
  
  matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
  matrixGENE <- rbind(matrixGENE, c(nGENE, length(iriGeneList$mixEnv)))
  if(length(gn_mixEnv) > 0) save(iriGeneList, file=paste0("Data/geneticsAS/retainedIntronByI_chr", chr, "_wt_p2.Rdata"))
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU

rownames(matrixGENE) <- paste0("chr", 1:5)
colnames(matrixGENE) <- c(paste0("Env", 1:4), "mixEnv")
matrixGENE
