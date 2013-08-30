#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 30-08-2013
# first written: 30-08-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#******************************************** this is the final version for analyzing results of testing retained introns ^_^ ********************************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for retained introns
riMatrix <- NULL
for(chr in 1:5){
  riMatrix <- rbind(riMatrix, read.table(paste0("Data/AS/retainedIntron_chr", chr, "_wt_p2.txt"), row.names=NULL))
}
#Bonferroni correction
nrow(riMatrix)# = 60981 introns were tested
-log10(0.05/nrow(riMatrix)/4)# = 6.69; 60981 introns were tested * 4 Env; => riThre=6.69
length(unique(riMatrix[,1]))# = 12995 expGenes have introns

#calculate the numbers of UN-sig introns and genes
riThre=round(-log10(0.05/82943/4), digits=2)# =6.82
matrixTU <- NULL #a matrix for numbers of introns that -log10(P) are lower than riThre in each env and across envs from chr1-chr5
matrixGENE <- NULL #a matrix for numbers of genes having at least one intron that -log10(P) are lower than riThre in each env and across envs from chr1-chr5
for(chr in 1:5){
  richr <- read.table(paste0("Data/AS/retainedIntron_chr", chr, "_wt_p2.txt"), row.names=NULL)
  riGeneList <- list()
  
  #FINDING UN-SIG ONES !!! 
  nTU <- NULL #number of introns that -log10(P) are lower than riThre in each env
  nGENE <- NULL #number of genes having one or more introns that -log10(P) are lower than riThre in each env
  for(env in 1:4){
    nTU <- c(nTU, length(which(richr[ ,env+2] < riThre)))
    riGeneList[[paste0("env", env)]] <- unique(richr[richr[ ,env+2] < riThre, 1]) #genes that have one or more introns that are spliced out in each env
    nGENE <- c(nGENE, length(riGeneList[[paste0("env", env)]]))
  }
  
  #cnt_mixEnv <- number of introns that -log10(P) are lower than riThre in ANY env
  cnt_mixEnv <- 0
  #gn_mixEnv <- genes having one or more introns that -log10(P) are lower than riThre in ANY env
  gn_mixEnv <- NULL
  for(e in 1:nrow(richr)){
    if(any(richr[e, 3:6] < riThre)){
      cnt_mixEnv <- cnt_mixEnv + 1
      gn_mixEnv <- c(gn_mixEnv, richr[e,1])
    }
  }
  riGeneList$mixEnv <- unique(gn_mixEnv)
  
  matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
  matrixGENE <- rbind(matrixGENE, c(nGENE, length(riGeneList$mixEnv)))
  if(length(gn_mixEnv) > 0) save(riGeneList, file=paste0("Data/AS/retainedIntron_chr", chr, "_wt_p2.Rdata"))
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU

rownames(matrixGENE) <- paste0("chr", 1:5)
colnames(matrixGENE) <- c(paste0("Env", 1:4), "mixEnv")
matrixGENE
