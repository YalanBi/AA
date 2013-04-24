#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 15-04-2013
# first written: 15-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]


#************************************************************ permutation part ************************************************************#
setwd("D:/Arabidopsis Arrays")
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]

#function for full model mapping
map.fast <- function(geno, pheno, menvironment){
  res <- NULL
  modelinfo <- summary(aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno) + as.factor(menvironment):as.numeric(geno)))
  res$qtl   <- unlist(lapply(modelinfo,"[",2,5),use.names=T)
  res$int   <- unlist(lapply(modelinfo,"[",3,5),use.names=T)
  res
}

#full model mapping with expressed genes
mapGenotypes <- function(geno, menvironment, permutation){
  st <- proc.time()
  cat("Start scanning for permutation", permutation, "\n")

  resQTL <- NULL
  resInt <- NULL
  
  for(m in 1:ncol(geno)){
    resList <- map.fast(geno[ ,m], pheno = expdata, menvironment)
    resQTL <- cbind(resQTL, -log10(resList$qtl))
    resInt <- cbind(resInt, -log10(resList$int))
  }
 
  et <- proc.time()
  cat("Done with full model mapping after:",(et-st)[3],"secs\n")
  return(c(max(resQTL), max(resInt)))
}

set.seed(1001)
#full model mapping, 5 chr
fileqtl <- "perms_qtls.txt"
fileint <- "perms_ints.txt"

cat("", file=fileqtl)
cat("", file=fileint)

chr <- 1
#the 1st col is bp!!!
expdata <- t(read.table(paste("Data/fullModeMapping/expGenes_chr", chr, "_1tu1probe.txt", sep=""),row.names=1, header=TRUE)[, 2:149])
expdata <- expdata[,1:100] # REMOVE BEFORE REAL !!!
  
for(permutation in 1:1000){
  newGeno <- geno[sample(1:148), ]
  scores <- mapGenotypes(geno = newGeno, menvironment, permutation)
  cat(paste0(permutation, "\t", scores[1],"\n"), file=fileqtl, append = TRUE)
  cat(paste0(permutation, "\t", scores[2],"\n"), file=fileint, append = TRUE)
}

