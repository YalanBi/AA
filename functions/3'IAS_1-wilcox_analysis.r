#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 09-07-2013
# first written: 28-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#****************************** this is the final version for analyzing results of testing GENETIC/INTERACTION regulated AS at 5/3 site ^_^ ******************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for INTERACTION regulated 3'AS
i3asMatrix <- NULL
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/3'AS_chr", chr, "_IAS_wt_less.txt"))){
    i3asMatrix <- rbind(i3asMatrix, read.table(paste0("Data/geneticsAS/3'AS_chr", chr, "_IAS_wt_less.txt"), row.names=NULL))
  } else cat("chr", chr, "NO test!\n")
}
#Bonferroni correction
-log10(0.05/nrow(i3asMatrix)/4)# = 1.90; 1 last exons were tested * 4 Env; => i3asThre=1.90
rm(i3asMatrix)

#calculate the numbers of sig exons and genes
i3asThre=1.9
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to i3asThre in each env and across envs from chr1-chr5
                 #NOTE: nExon=nGene(one first exon in one gene)
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/3'AS_chr", chr, "_IAS_wt_less.txt"))){
    i3aschr <- read.table(paste0("Data/geneticsAS/3'AS_chr", chr, "_IAS_wt_less.txt"), row.names=NULL)
    i3asGeneList <- list()
    resmatrix <- NULL
    
    for(e in 1:nrow(i3aschr)){
      res <- NULL
      for(env in 1:4){
        if((i3aschr[e,env*2+2] >= i3asThre && i3aschr[e,env*2+3] < i3asThre) || (i3aschr[e,env*2+2] < i3asThre && i3aschr[e,env*2+3] >= i3asThre)){
          res <- c(res, 1)
        } else res <- c(res, 0)
      }
      resmatrix <- rbind(resmatrix, res)
    }
    
    nTU <- NULL #number of exons that -log10(P) are higher than or equal to i3asThre in one gt and lower in the other gt, each env sep
    for(env in 1:4){
      nTU <- c(nTU, sum(resmatrix[ ,env]))
      i3asGeneList[[paste0("env", env)]] <- unique(i3aschr[resmatrix[ ,env] > 0, 1]) #genes that have one or more exons that are spliced out in each env
    }
    
    #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to i3asThre in ANY env
    cnt_mixEnv <- 0
    #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to i3asThre in ANY env
    gn_mixEnv <- NULL
    for(e in 1:nrow(i3aschr)){
      if(sum(resmatrix[e, ]) > 0){
        cnt_mixEnv <- cnt_mixEnv + 1
        gn_mixEnv <- c(gn_mixEnv, i3aschr[e,1])
      }
    }
    i3asGeneList$mixEnv <- unique(gn_mixEnv)
    
    matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
    if(length(gn_mixEnv) > 0) save(i3asGeneList, file = paste0("Data/geneticsAS/3'AS_chr", chr, "_IAS_gnList_wt.Rdata"))
    else cat("chr", chr, "NO genes, no file saved!\n")
  } else{
    cat("chr", chr, "NO test!\n")
    matrixTU <- rbind(matrixTU, c(0,0,0,0,0))
  }
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU
