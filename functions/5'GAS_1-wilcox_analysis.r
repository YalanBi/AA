#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 09-07-2013
# first written: 28-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#****************************** this is the final version for analyzing results of testing GENETIC/INTERACTION regulated AS at 5/3 site ^_^ ******************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for GENETIC regulated 5'AS
g5asMatrix <- NULL
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/5'AS_chr", chr, "_GAS_wt_less.txt"))){
    g5asMatrix <- rbind(g5asMatrix, read.table(paste0("Data/geneticsAS/5'AS_chr", chr, "_GAS_wt_less.txt"), row.names=NULL))
  } else cat("chr", chr, "NO test!\n")
}
#Bonferroni correction
-log10(0.05/nrow(g5asMatrix)/4)# = 2.68; 6 first exons were tested * 4 Env; => g5asThre=2.68
rm(g5asMatrix)

#calculate the numbers of sig exons and genes
g5asThre=2.68
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to g5asThre in each env and across envs from chr1-chr5
                 #NOTE: nExon=nGene(one first exon in one gene)
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/5'AS_chr", chr, "_GAS_wt_less.txt"))){
    g5aschr <- read.table(paste0("Data/geneticsAS/5'AS_chr", chr, "_GAS_wt_less.txt"), row.names=NULL)
    g5asGeneList <- list()
    resmatrix <- NULL
    
    for(e in 1:nrow(g5aschr)){
      res <- NULL
      for(env in 1:4){
        if((g5aschr[e,env*2+2] >= g5asThre && g5aschr[e,env*2+3] < g5asThre) || (g5aschr[e,env*2+2] < g5asThre && g5aschr[e,env*2+3] >= g5asThre)){
          res <- c(res, 1)
        } else res <- c(res, 0)
      }
      resmatrix <- rbind(resmatrix, res)
    }
    
    nTU <- NULL #number of exons that -log10(P) are higher than or equal to g5asThre in one gt and lower in the other gt, each env sep
    for(env in 1:4){
      nTU <- c(nTU, sum(resmatrix[ ,env]))
      g5asGeneList[[paste0("env", env)]] <- unique(g5aschr[resmatrix[ ,env] > 0, 1]) #genes that have one or more exons that are spliced out in each env
    }
    
    #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to g5asThre in ANY env
    cnt_mixEnv <- 0
    #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to g5asThre in ANY env
    gn_mixEnv <- NULL
    for(e in 1:nrow(g5aschr)){
      if(sum(resmatrix[e, ]) > 0){
        cnt_mixEnv <- cnt_mixEnv + 1
        gn_mixEnv <- c(gn_mixEnv, g5aschr[e,1])
      }
    }
    g5asGeneList$mixEnv <- unique(gn_mixEnv)
    
    matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
    if(length(gn_mixEnv) > 0) save(g5asGeneList, file = paste0("Data/geneticsAS/5'AS_chr", chr, "_GAS_gnList_wt.Rdata"))
    else cat("chr", chr, "NO genes, no file saved!\n")
  } else{
    cat("chr", chr, "NO test!\n")
    matrixTU <- rbind(matrixTU, c(0,0,0,0,0))
  }
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU
