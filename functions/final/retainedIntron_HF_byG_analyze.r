#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 03-09-2013
# first written: 03-09-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#*********************************** this is the final version for analyzing results of testing GENETIC regulated retained introns ^_^ ***********************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for GENETIC regulated retained introns
griMatrix <- NULL
for(chr in 1:5){
  griMatrix <- rbind(griMatrix, read.table(paste0("Data/geneticsAS/retainedIntronByG_chr", chr, "_wt_p2.txt"), row.names=NULL))
}
#Bonferroni correction
nrow(griMatrix)# = 60981 introns were tested 
-log10(0.05/nrow(griMatrix)/8)# 60981 introns were tested * 8 Env; => griThre=6.99
length(unique(griMatrix[,1]))# = 12995 genes were tested

#calculate the numbers of NOT sig introns and genes
griThre=round(-log10(0.05/82943/8), digits=2)# =7.12
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to griThre in each env and across envs from chr1-chr5
matrixGENE <- NULL #a matrix for numbers of genes having at least one exon that -log10(P) are higher than or equal to griThre in each env and across envs from chr1-chr5
for(chr in 1:5){
  grichr <- read.table(paste0("Data/geneticsAS/retainedIntronByG_chr", chr, "_wt_p2.txt"), row.names=NULL)
  griGeneList <- list()
  resmatrix <- NULL
  
  for(e in 1:nrow(grichr)){
    res <- NULL
    for(env in 1:4){
      if((grichr[e,env*2+2] >= griThre && grichr[e,env*2+3] < griThre) || (grichr[e,env*2+2] < griThre && grichr[e,env*2+3] >= griThre)){
        res <- c(res, 1)
      }else res <- c(res, 0)
    }
    resmatrix <- rbind(resmatrix, res)
  }
  
  nTU <- NULL #number of exons that -log10(P) are higher than or equal to griThre in one gt and lower in the other gt, each env sep
  nGENE <- NULL #number of genes having one or more exons that -log10(P) are higher than or equal to griThre in each env
  for(env in 1:4){
    nTU <- c(nTU, sum(resmatrix[ ,env]))
    griGeneList[[paste0("env", env)]] <- unique(grichr[resmatrix[ ,env] > 0, 1]) #genes that have one or more exons that are spliced out in each env
    nGENE <- c(nGENE, length(griGeneList[[paste0("env", env)]]))
  }
  
  #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to griThre in ANY env
  cnt_mixEnv <- 0
  #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to griThre in ANY env
  gn_mixEnv <- NULL
  for(e in 1:nrow(grichr)){
    if(sum(resmatrix[e, ]) > 0){
      cnt_mixEnv <- cnt_mixEnv + 1
      gn_mixEnv <- c(gn_mixEnv, grichr[e,1])
    }
  }
  griGeneList$mixEnv <- unique(gn_mixEnv)
  
  matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
  matrixGENE <- rbind(matrixGENE, c(nGENE, length(griGeneList$mixEnv)))
  if(length(gn_mixEnv) > 0) save(griGeneList, file=paste0("Data/geneticsAS/retainedIntronByG_chr", chr, "_wt_p2.Rdata"))
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU

rownames(matrixGENE) <- paste0("chr", 1:5)
colnames(matrixGENE) <- c(paste0("Env", 1:4), "mixEnv")
matrixGENE
