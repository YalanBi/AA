#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 20-08-2013
# first written: 28-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#************************************ this is the final version for analyzing results of testing GENETIC regulated skipping exons ^_^ ************************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for GENETIC regulated skipping exons
gseMatrix <- NULL
for(chr in 1:5){
  gseMatrix <- rbind(gseMatrix, read.table(paste0("Data/geneticsAS/skippingExonByG_chr", chr, "_wt_p3.txt"), row.names=NULL))
}
#Bonferroni correction
nrow(gseMatrix)# = 98677 exons were tested
-log10(0.05/nrow(gseMatrix)/8)# = 7.20; 98677 exons were tested * 4 Env * 2 genotypes; => gseThre=7.20
length(unique(gseMatrix[,1]))# = 16448 genes were tested

#calculate the numbers of sig exons and genes
gseThre=round(-log10(0.05/nrow(gseMatrix)/8), digits=2)# =7.20
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to gseThre in each env and across envs from chr1-chr5
matrixGENE <- NULL #a matrix for numbers of genes having at least one exon that -log10(P) are higher than or equal to gseThre in each env and across envs from chr1-chr5
for(chr in 1:5){
  gsechr <- read.table(paste0("Data/geneticsAS/skippingExonByG_chr", chr, "_wt_p3.txt"), row.names=NULL)
  gseGeneList <- list()
  resmatrix <- NULL
  
  for(e in 1:nrow(gsechr)){
    res <- NULL
    for(env in 1:4){
      if((gsechr[e,env*2+2] >= gseThre && gsechr[e,env*2+3] < gseThre) || (gsechr[e,env*2+2] < gseThre && gsechr[e,env*2+3] >= gseThre)){
        res <- c(res, 1)
      } else res <- c(res, 0)
    }
    resmatrix <- rbind(resmatrix, res)
  }
  
  nTU <- NULL #number of exons that -log10(P) are higher than or equal to gseThre in one gt and lower in the other gt, each env sep
  nGENE <- NULL #number of genes having one or more exons that -log10(P) are higher than or equal to gseThre in each env
  for(env in 1:4){
    nTU <- c(nTU, sum(resmatrix[ ,env]))
    gseGeneList[[paste0("env", env)]] <- unique(gsechr[resmatrix[ ,env] > 0, 1]) #genes that have one or more exons that are spliced out in each env
    nGENE <- c(nGENE, length(gseGeneList[[paste0("env", env)]]))
  }
  
  #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to gseThre in ANY env
  cnt_mixEnv <- 0
  #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to gseThre in ANY env
  gn_mixEnv <- NULL
  for(e in 1:nrow(gsechr)){
    if(sum(resmatrix[e, ]) > 0){
      cnt_mixEnv <- cnt_mixEnv + 1
      gn_mixEnv <- c(gn_mixEnv, gsechr[e,1])
    }
  }
  gseGeneList$mixEnv <- unique(gn_mixEnv)
  
  matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
  matrixGENE <- rbind(matrixGENE, c(nGENE, length(gseGeneList$mixEnv)))
  if(length(gn_mixEnv) > 0) save(gseGeneList, file=paste0("Data/geneticsAS/skippingExonByG_chr", chr, "_wt_p3.Rdata"))
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU

rownames(matrixGENE) <- paste0("chr", 1:5)
colnames(matrixGENE) <- c(paste0("Env", 1:4), "mixEnv")
matrixGENE
