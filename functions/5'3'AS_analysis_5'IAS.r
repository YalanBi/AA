#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 01-07-2013
# first written: 28-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#****************************** this is the final version for analyzing results of testing GENETIC/INTERACTION regulated AS at 5/3 site ^_^ ******************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for INTERACTION regulated 5'AS
i5asMatrix <- NULL
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/5'AS_chr", chr, "_IAS_wt_less.txt"))){
    i5asMatrix <- rbind(i5asMatrix, read.table(paste0("Data/geneticsAS/5'AS_chr", chr, "_IAS_wt_less.txt"), row.names=NULL))
  } else cat("chr", chr, "NO test!\n")
}
#Bonferroni correction
-log10(0.05/nrow(i5asMatrix)/4)# = 2.51; 4 first exons were tested * 4 Env; => i5asThre=2.51
rm(i5asMatrix)

#calculate the numbers of sig exons and genes
i5asThre=2.51
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to i5asThre in each env and across envs from chr1-chr5
                 #NOTE: nExon=nGene(one first exon in one gene)
for(chr in 1:5){
  if(file.exists(paste0("Data/geneticsAS/5'AS_chr", chr, "_IAS_wt_less.txt"))){
    i5aschr <- read.table(paste0("Data/geneticsAS/5'AS_chr", chr, "_IAS_wt_less.txt"), row.names=NULL)
    i5asGeneList <- list()
    resmatrix <- NULL
    
    for(e in 1:nrow(i5aschr)){
      res <- NULL
      for(env in 1:4){
        if((i5aschr[e,env*2+2] >= i5asThre && i5aschr[e,env*2+3] < i5asThre) || (i5aschr[e,env*2+2] < i5asThre && i5aschr[e,env*2+3] >= i5asThre)){
          res <- c(res, 1)
        } else res <- c(res, 0)
      }
      resmatrix <- rbind(resmatrix, res)
    }
    
    nTU <- NULL #number of exons that -log10(P) are higher than or equal to i5asThre in one gt and lower in the other gt, each env sep
    for(env in 1:4){
      nTU <- c(nTU, sum(resmatrix[ ,env]))
      i5asGeneList[[paste0("env", env)]] <- unique(i5aschr[resmatrix[ ,env] > 0, 1]) #genes that have one or more exons that are spliced out in each env
    }
    
    #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to i5asThre in ANY env
    cnt_mixEnv <- 0
    #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to i5asThre in ANY env
    gn_mixEnv <- NULL
    for(e in 1:nrow(i5aschr)){
      if(sum(resmatrix[e, ]) > 0){
        cnt_mixEnv <- cnt_mixEnv + 1
        gn_mixEnv <- c(gn_mixEnv, i5aschr[e,1])
      }
    }
    i5asGeneList$mixEnv <- unique(gn_mixEnv)
    
    matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
    if(length(gn_mixEnv) > 0) save(i5asGeneList, file = paste0("Data/geneticsAS/5'AS_chr", chr, "_IAS_genenameList.Rdata"))
    else cat("chr", chr, "NO genes, no file saved!\n")
  } else{
    cat("chr", chr, "NO test!\n")
    matrixTU <- rbind(matrixTU, c(0,0,0,0,0))
  }
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU
