#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 27-08-2013
# first written: 22-08-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#************************************************ this is the final version for analyzing results of testing 5'3'AS ! ^_^ ************************************************#
#**************************************************************** testing algorithm: Wilcox.test! ****************************************************************#

setwd("D:/Arabidopsis Arrays")

#calculate the threshold for 3'AS
asMatrix <- NULL
for(chr in 1:5){
  asMatrix <- rbind(asMatrix, read.table(paste0("Data/AS/splicing3'site_chr", chr, "_wt_p2.txt"), row.names=NULL))
}
#Bonferroni correction
nrow(asMatrix)# = 16648 exons were tested in 3'AS
length(unique(asMatrix[ ,1]))# = 16648 genes were tested in 3'AS
-log10(0.05/nrow(asMatrix)/4)# 16648 exons were tested * 4 Env; => asThre=6.12

#calculate the numbers of sig exons and genes
asThre=round(-log10(0.05/30528/4), digits=2)# =6.39
matrixTU <- NULL #a matrix for numbers of exons that -log10(P) are higher than or equal to asThre in each env and across envs from chr1-chr5,
                 #NOTE: nExon=nGene(one first exon in one gene)
for(chr in 1:5){
  aschr <- read.table(paste0("Data/AS/splicing3'site_chr", chr, "_wt_p2.txt"), row.names=NULL)
  asGeneList <- list()
  
  nTU <- NULL #number of exons that -log10(P) are higher than or equal to asThre in each env
  for(env in 1:4){
    nTU <- c(nTU, length(which(aschr[ ,env+2] >= asThre))) #now, the matrix has no rowname, so env+2
    asGeneList[[paste0("env", env)]] <- aschr[aschr[ ,env+2] >= asThre, 1] #genes that its first exon is spliced out in each env separately
  }
  
  #cnt_mixEnv <- number of exons that -log10(P) are higher than or equal to asThre in ANY env
  cnt_mixEnv <- 0
  #gn_mixEnv <- genes having one or more exons that -log10(P) are higher than or equal to asThre in ANY env
  gn_mixEnv <- NULL
  for(e in 1:nrow(aschr)){
    if(any(aschr[e, 3:6] >= asThre)){
      cnt_mixEnv <- cnt_mixEnv + 1
      gn_mixEnv <- c(gn_mixEnv, aschr[e, 1])
    }
  }
  asGeneList$mixEnv <- gn_mixEnv
  
  matrixTU <- rbind(matrixTU, c(nTU, cnt_mixEnv))
  if(length(gn_mixEnv) > 0) save(asGeneList, file = paste0("Data/AS/splicing3'site_chr", chr, "_wt_p2.Rdata"))
  else cat("chr", chr, "NO genes, no file saved!\n")
}
rownames(matrixTU) <- paste0("chr", 1:5)
colnames(matrixTU) <- c(paste0("Env", 1:4), "mixEnv")
matrixTU
