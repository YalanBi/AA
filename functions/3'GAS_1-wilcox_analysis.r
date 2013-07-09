#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 01-07-2013
# first written: 28-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#****************************** this is the final version for analyzing results of testing GENETIC/INTERACTION regulated AS at 5/3 site ^_^ ******************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for GENETIC regulated 3'AS
g3asMatrix <- NULL
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/3'AS_chr", chr, "_GAS_wt_less.txt"))){
    g3asMatrix <- rbind(g3asMatrix, read.table(paste0("Data/geneticsAS/3'AS_chr", chr, "_GAS_wt_less.txt"), row.names=NULL))
  } else cat("chr", chr, "NO test!\n")
}
#Bonferroni correction
-log10(0.05/nrow(g3asMatrix)/4)# = 2.98; 12 last exons were tested * 4 Env; => g3asThre=2.98
rm(g3asMatrix)

#calculate the numbers of sig exons and genes
g3asThre=2.98
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to g3asThre in each env and across envs from chr1-chr5
                 #NOTE: nExon=nGene(one first exon in one gene)
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/3'AS_chr", chr, "_GAS_wt_less.txt"))){
    g3aschr <- read.table(paste0("Data/geneticsAS/3'AS_chr", chr, "_GAS_wt_less.txt"), row.names=NULL)
    g3asGeneList <- list()
    resmatrix <- NULL
    
    for(e in 1:nrow(g3aschr)){
      res <- NULL
      for(env in 1:4){
        if((g3aschr[e,env*2+2] >= g3asThre && g3aschr[e,env*2+3] < g3asThre) || (g3aschr[e,env*2+2] < g3asThre && g3aschr[e,env*2+3] >= g3asThre)){
          res <- c(res, 1)
        } else res <- c(res, 0)
      }
      resmatrix <- rbind(resmatrix, res)
    }
    
    nTU <- NULL #number of exons that -log10(P) are higher than or equal to g3asThre in one gt and lower in the other gt, each env sep
    for(env in 1:4){
      nTU <- c(nTU, sum(resmatrix[ ,env]))
      g3asGeneList[[paste0("env", env)]] <- unique(g3aschr[resmatrix[ ,env] > 0, 1]) #genes that have one or more exons that are spliced out in each env
    }
    
    #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to g3asThre in ANY env
    cnt_mixEnv <- 0
    #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to g3asThre in ANY env
    gn_mixEnv <- NULL
    for(e in 1:nrow(g3aschr)){
      if(sum(resmatrix[e, ]) > 0){
        cnt_mixEnv <- cnt_mixEnv + 1
        gn_mixEnv <- c(gn_mixEnv, g3aschr[e,1])
      }
    }
    g3asGeneList$mixEnv <- unique(gn_mixEnv)
    
    matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
    if(length(gn_mixEnv) > 0) save(g3asGeneList, file = paste0("Data/geneticsAS/3'AS_chr", chr, "_GAS_gnList_wt.Rdata"))
    else cat("chr", chr, "NO genes, no file saved!\n")
  } else{
    cat("chr", chr, "NO test!\n")
    matrixTU <- rbind(matrixTU, c(0,0,0,0,0))
  }
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU
