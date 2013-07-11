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

#calculate the threshold for INTERACTION regulated retained introns
iriMatrix <- NULL
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/IRI_chr", chr, "_", whichTest, "_less.txt"))){
    iriMatrix <- rbind(iriMatrix, read.table(paste0("Data/geneticsAS/IRI_chr", chr, "_", whichTest, "_less.txt"), row.names=NULL))
  } else cat("chr", chr, "NO test!\n")
}
#Bonferroni correction
nrow(iriMatrix)# = 1401 exons were tested in 5'GAS; 1525 exons were tested in 3'GAS
-log10(0.05/nrow(iriMatrix)/4)# 1401/1525 exons were tested * 4 Env; => iriThre=5.05/5.09
length(unique(iriMatrix[,1]))# = 1401 genes were tested in 5'GAS; 1525 genes were tested in 3'GAS
#STOP!



#************************************************************************** NO test on any chr! **************************************************************************#
#Bonferroni correction
-log10(0.05/nrow(iriMatrix)/4)# = 2.51; 4 introns were tested * 4 Env; => iriThre=2.51
rm(iriMatrix)

#calculate the numbers of sig exons and genes
iriThre=2.51
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to iriThre in each env and across envs from chr1-chr5
matrixGENE <- NULL #a matrix for numbers of genes having at least one exon that -log10(P) are higher than or equal to iriThre in each env and across envs from chr1-chr5
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/RI_chr", chr, "_GAS_wt_less.txt"))){
    irichr <- read.table(paste0("Data/geneticsAS/RI_chr", chr, "_GAS_wt_less.txt"), row.names=NULL)
    iriGeneList <- list()
    resmatrix <- NULL
    
    for(e in 1:nrow(irichr)){
      res <- NULL
      for(env in 1:4){
        if((irichr[e,env*2+2] >= iriThre && irichr[e,env*2+3] < iriThre) || (irichr[e,env*2+2] < iriThre && irichr[e,env*2+3] >= iriThre)){
          res <- c(res, 1)
        } else res <- c(res, 0)
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
    save(iriGeneList, file = paste0("Data/geneticsAS/RI_chr", chr, "_GAS_genenameList.Rdata"))
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
