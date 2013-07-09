#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 09-07-2013
# first written: 21-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#********************************************* this is the final version for analyzing results of testing AS at 3' site! ^_^ *********************************************#
#*********************************************************************** testing algorithm: ANOVA! ***********************************************************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for 3'AS
as3Matrix <- NULL
for(chr in 1:5){
  as3Matrix <- rbind(as3Matrix, read.table(paste0("Data/AS/3AS_chr", chr, "_ANOVA_less.txt"), row.names=NULL))
}
#Bonferroni correction
-log10(0.05/nrow(as3Matrix)/4)# = 4.63; 534 exons were tested * 4 Env; => as3Thre=4.63
rm(as3Matrix)

#calculate the numbers of sig exons and genes
as3Thre=4.63
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to as3Thre in each env and across envs from chr1-chr5,
                 #NOTE: nExon=nGene(one first exon in one gene)
for(chr in 1:5){
  as3chr <- read.table(paste0("Data/AS/3AS_chr", chr, "_ANOVA_less.txt"), row.names=NULL)
  as3GeneList <- list()
  
  nTU <- NULL #number of exons that -log10(P) are higher than or equal to as3Thre in each env
  for(env in 1:4){
    nTU <- c(nTU, length(which(as3chr[ ,env+2] >= as3Thre))) #now, the matrix has no rowname, so env+2
    as3GeneList[[paste0("env", env)]] <- as3chr[as3chr[ ,env+2] >= as3Thre, 1] #genes that its first exon is spliced out in each env separately
  }
  
  #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to as3Thre in ANY env
  cnt_mixEnv <- 0
  #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to as3Thre in ANY env
  gn_mixEnv <- NULL
  for(e in 1:nrow(as3chr)){
    if(any(as3chr[e, 3:6] >= as3Thre)){
      cnt_mixEnv <- cnt_mixEnv + 1
      gn_mixEnv <- c(gn_mixEnv, as3chr[e, 1])
    }
  }
  as3GeneList$mixEnv <- gn_mixEnv
  
  matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
  save(as3GeneList, file = paste0("Data/AS/3AS_chr", chr, "_gnList_ANOVA.Rdata"))
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU
