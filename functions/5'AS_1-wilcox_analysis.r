#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 09-07-2013
# first written: 21-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#********************************************* this is the final version for analyzing results of testing AS at 5' site! ^_^ *********************************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for 5'AS
as5Matrix <- NULL
for(chr in 1:5){
  as5Matrix <- rbind(as5Matrix, read.table(paste0("Data/AS/5AS_chr", chr, "_wt_less.txt"), row.names=NULL))
}
#Bonferroni correction
-log10(0.05/nrow(as5Matrix)/4)# = 4.54; 429 exons were tested * 4 Env; => as5Thre=4.54
rm(as5Matrix)

#calculate the numbers of sig exons and genes
as5Thre=4.54
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to as5Thre in each env and across envs from chr1-chr5,
                 #NOTE: nExon=nGene(one first exon in one gene)
for(chr in 1:5){
  as5chr <- read.table(paste0("Data/AS/5AS_chr", chr, "_wt_less.txt"), row.names=NULL)
  as5GeneList <- list()
  
  nTU <- NULL #number of exons that -log10(P) are higher than or equal to as5Thre in each env
  for(env in 1:4){
    nTU <- c(nTU, length(which(as5chr[ ,env+2] >= as5Thre))) #now, the matrix has no rowname, so env+2
    as5GeneList[[paste0("env", env)]] <- as5chr[as5chr[ ,env+2] >= as5Thre, 1] #genes that its first exon is spliced out in each env separately
  }
  
  #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to as5Thre in ANY env
  cnt_mixEnv <- 0
  #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to as5Thre in ANY env
  gn_mixEnv <- NULL
  for(e in 1:nrow(as5chr)){
    if(any(as5chr[e, 3:6] >= as5Thre)){
      cnt_mixEnv <- cnt_mixEnv + 1
      gn_mixEnv <- c(gn_mixEnv, as5chr[e, 1])
    }
  }
  as5GeneList$mixEnv <- gn_mixEnv
  
  matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
  save(as5GeneList, file = paste0("Data/AS/5AS_chr", chr, "_gnList_wt.Rdata"))
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU
