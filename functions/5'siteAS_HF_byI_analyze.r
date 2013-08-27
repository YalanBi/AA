#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 27-08-2013
# first written: 27-08-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#************************************ this is the final version for analyzing results of testing interaction regulated AS 5'site ^_^ *************************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for interaction regulated 5'site alternative splicing
iasMatrix <- NULL
for(chr in 1:5){
  iasMatrix <- rbind(iasMatrix, read.table(paste0("Data/geneticsAS/splicing5'siteByI_chr", chr, "_wt_p2.txt"), row.names=NULL))
}
#Bonferroni correction
nrow(iasMatrix)# = 16448 exons were tested
-log10(0.05/nrow(iasMatrix)/8)# = 6.42; 16448 exons were tested * 4 Env * 2 genotypes
length(unique(iasMatrix[,1]))# = 16448 genes were tested

#calculate the numbers of sig exons and genes
iasThre=round(-log10(0.05/30528/8), digits=2)# = 6.69
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to iasThre in each env and across envs from chr1-chr5
                 #NOTE: nExon=nGene(one first exon in one gene)
for(chr in 1:5){
  iaschr <- read.table(paste0("Data/geneticsAS/splicing5'siteByI_chr", chr, "_wt_p2.txt"), row.names=NULL)
  iasGeneList <- list()
  resmatrix <- NULL
  
  for(e in 1:nrow(iaschr)){
    res <- NULL
    for(env in 1:4){
      if((iaschr[e,env*2+2] >= iasThre && iaschr[e,env*2+3] < iasThre) || (iaschr[e,env*2+2] < iasThre && iaschr[e,env*2+3] >= iasThre)){
        res <- c(res, 1)
      } else res <- c(res, 0)
    }
    resmatrix <- rbind(resmatrix, res)
  }
  
  nTU <- NULL #number of exons that -log10(P) are higher than or equal to iasThre in one gt and lower in the other gt, each env sep
  for(env in 1:4){
    nTU <- c(nTU, sum(resmatrix[ ,env]))
    iasGeneList[[paste0("env", env)]] <- unique(iaschr[resmatrix[ ,env] > 0, 1]) #genes that have one or more exons that are spliced out in each env
  }
  
  #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to iasThre in ANY env
  cnt_mixEnv <- 0
  #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to iasThre in ANY env
  gn_mixEnv <- NULL
  for(e in 1:nrow(iaschr)){
    if(sum(resmatrix[e, ]) > 0){
      cnt_mixEnv <- cnt_mixEnv + 1
      gn_mixEnv <- c(gn_mixEnv, iaschr[e,1])
    }
  }
  iasGeneList$mixEnv <- unique(gn_mixEnv)
  
  matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
  if(length(gn_mixEnv) > 0) save(iasGeneList, file=paste0("Data/geneticsAS/splicing5'siteByI_chr", chr, "_wt_p2.Rdata"))
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU
