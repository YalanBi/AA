#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 12-07-2013
# first written: 28-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#****************************** this is the final version for analyzing results of testing GENETIC/INTERACTION regulated AS at 5/3 site ^_^ ******************************#
#**************************************************************** testing algorithm: Wilcox.test / ANOVA! ****************************************************************#

setwd("D:/Arabidopsis Arrays")

goal="5'IAS"# "3'IAS"
whichTest="wt"# "ANOVA"

#calculate the threshold for INTERACTION regulated 5'3'AS
iasMatrix <- NULL
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/", goal, "_chr", chr, "_", whichTest, "_less.txt"))){
    iasMatrix <- rbind(iasMatrix, read.table(paste0("Data/geneticsAS/", goal, "_chr", chr, "_", whichTest, "_less.txt"), row.names=NULL))
  } else cat("chr", chr, "NO test!\n")
}
#Bonferroni correction
nrow(iasMatrix)# = 4 exons were tested in 5'IAS; 6 exons were tested in 3'IAS
-log10(0.05/nrow(iasMatrix)/4)# 4/6 exons were tested * 4 Env; => iasThre=2.51/2.68
length(unique(iasMatrix[,1]))# = 4 genes were tested in 5'IAS; 6 genes were tested in 3'IAS

#calculate the numbers of sig exons and genes
iasThre=round(-log10(0.05/nrow(iasMatrix)/4), digits=2)# =2.51/2.68
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to iasThre in each env and across envs from chr1-chr5
                 #NOTE: nExon=nGene(one first exon in one gene)
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/", goal, "_chr", chr, "_", whichTest, "_less.txt"))){
    iaschr <- read.table(paste0("Data/geneticsAS/", goal, "_chr", chr, "_", whichTest, "_less.txt"), row.names=NULL)
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
    if(length(gn_mixEnv) > 0) save(iasGeneList, file = paste0("Data/geneticsAS/", goal, "_chr", chr, "_", whichTest, ".Rdata"))
    else cat("chr", chr, "NO genes, no file saved!\n")
  } else{
    cat("chr", chr, "NO test!\n")
    matrixTU <- rbind(matrixTU, c(0,0,0,0,0))
  }
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU
