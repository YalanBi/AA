#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 30-08-2013
# first written: 30-08-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#************************************** this is the final version for analyzing results of testing genetics regulated AS 3'site ^_^ **************************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for genetics regulated 3'site alternative splicing
gasMatrix <- NULL
for(chr in 1:5){
  gasMatrix <- rbind(gasMatrix, read.table(paste0("Data/geneticsAS/splicing3'siteByG_chr", chr, "_wt_p2.txt"), row.names=NULL))
}
#Bonferroni correction
nrow(gasMatrix)# = 16448 exons were tested
-log10(0.05/nrow(gasMatrix)/8)# = 6.42; 16448 exons were tested * 4 Env * 2 genotypes
length(unique(gasMatrix[,1]))# = 16448 genes were tested

#calculate the numbers of sig exons and genes
gasThre=round(-log10(0.05/30528/8), digits=2)# = 6.69
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to gasThre in each env and across envs from chr1-chr5
                 #NOTE: nExon=nGene(one first exon in one gene)
for(chr in 1:5){
  gaschr <- read.table(paste0("Data/geneticsAS/splicing3'siteByG_chr", chr, "_wt_p2.txt"), row.names=NULL)
  gasGeneList <- list()
  resmatrix <- NULL
  
  for(e in 1:nrow(gaschr)){
    res <- NULL
    for(env in 1:4){
      if((gaschr[e,env*2+2] >= gasThre && gaschr[e,env*2+3] < gasThre) || (gaschr[e,env*2+2] < gasThre && gaschr[e,env*2+3] >= gasThre)){
        res <- c(res, 1)
      } else res <- c(res, 0)
    }
    resmatrix <- rbind(resmatrix, res)
  }
  
  nTU <- NULL #number of exons that -log10(P) are higher than or equal to gasThre in one gt and lower in the other gt, each env sep
  for(env in 1:4){
    nTU <- c(nTU, sum(resmatrix[ ,env]))
    gasGeneList[[paste0("env", env)]] <- unique(gaschr[resmatrix[ ,env] > 0, 1]) #genes that have one or more exons that are spliced out in each env
  }
  
  #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to gasThre in ANY env
  cnt_mixEnv <- 0
  #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to gasThre in ANY env
  gn_mixEnv <- NULL
  for(e in 1:nrow(gaschr)){
    if(sum(resmatrix[e, ]) > 0){
      cnt_mixEnv <- cnt_mixEnv + 1
      gn_mixEnv <- c(gn_mixEnv, gaschr[e,1])
    }
  }
  gasGeneList$mixEnv <- unique(gn_mixEnv)
  
  matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
  if(length(gn_mixEnv) > 0) save(gasGeneList, file=paste0("Data/geneticsAS/splicing3'siteByG_chr", chr, "_wt_p2.Rdata"))
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU
