#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 11-07-2013
# first written: 28-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#***************************** this is the final version for analyzing results of testing GENETIC/INTERACTION regulated retained introns ^_^ *****************************#
#**************************************************************** testing algorithm: Wilcox.test / ANOVA! ****************************************************************#

setwd("D:/Arabidopsis Arrays")

whichTest="wt"# "ANOVA"

#calculate the threshold for GENETIC regulated retained introns
griMatrix <- NULL
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/GRI_chr", chr, "_", whichTest, "_less.txt"))){
    griMatrix <- rbind(griMatrix, read.table(paste0("Data/geneticsAS/GRI_chr", chr, "_", whichTest, "_less.txt"), row.names=NULL))
  } else cat("chr", chr, "NO test!\n")
}
#Bonferroni correction
nrow(griMatrix)# = 1401 exons were tested in 5'GAS; 1525 exons were tested in 3'GAS
-log10(0.05/nrow(griMatrix)/4)# 1401/1525 exons were tested * 4 Env; => griThre=5.05/5.09
length(unique(griMatrix[,1]))# = 1401 genes were tested in 5'GAS; 1525 genes were tested in 3'GAS

#calculate the numbers of sig exons and genes
griThre=round(-log10(0.05/nrow(griMatrix)/4), digits=2)# =5.05/5.09
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to griThre in each env and across envs from chr1-chr5
matrixGENE <- NULL #a matrix for numbers of genes having at least one exon that -log10(P) are higher than or equal to griThre in each env and across envs from chr1-chr5
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/GRI_chr", chr, "_", whichTest, "_less.txt"))){
    grichr <- read.table(paste0("Data/geneticsAS/GRI_chr", chr, "_", whichTest, "_less.txt"), row.names=NULL)
    griGeneList <- list()
    resmatrix <- NULL
    
    for(e in 1:nrow(grichr)){
      res <- NULL
      for(env in 1:4){
        if((grichr[e,env*2+2] >= griThre && grichr[e,env*2+3] < griThre) || (grichr[e,env*2+2] < griThre && grichr[e,env*2+3] >= griThre)){
          res <- c(res, 1)
        } else res <- c(res, 0)
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
    save(griGeneList, file = paste0("Data/geneticsAS/GRI_chr", chr, "_gnList_", whichTest, ".Rdata"))
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
