#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 22-03-2013
# first written: 22-03-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("D:/Arabidopsis Arrays")

geno <- read.table("refined map/genotypes.txt", row.names=1)
env <- read.table("Data/ann_env.txt",sep="\t")[,2]

map.fast <- function(geno, pheno, menvironment){
  res <- NULL
  modelinfo <- summary(aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno) + as.factor(menvironment):as.numeric(geno)))
  res$env   <- unlist(lapply(modelinfo,"[",1,5), use.names=T)
  res$qtl   <- unlist(lapply(modelinfo,"[",2,5), use.names=T)
  res$int   <- unlist(lapply(modelinfo,"[",3,5), use.names=T)
  return(res)
}

mapGenesByTU <- function(filename = "genesByTU_chr1_norm_hf_cor.txt", geno, env){
  location <- "D:/Arabidopsis Arrays/Data/"
  st <- proc.time()[3]
  new_exp <- t(read.table(paste("Data/", filename, sep=""), row.names=1, header=TRUE))
  resEnv <- NULL
  resQTL <- NULL
  resInt <- NULL
  
  resEnv <- apply(geno, 2, function(x, pheno, env){
    -log10(map.fast(x, pheno, env)$env)
  }, pheno=new_exp, env=env)
  rownames(resEnv) <- colnames(new_exp)
  colnames(resEnv) <- colnames(geno)
  write.table(resEnv, file=paste(location, gsub(".txt", "_Env.txt", filename), sep=""))
  
  resQTL <- apply(geno, 2, function(x, pheno, env){
    -log10(map.fast(x, pheno, env)$qtl)
  }, pheno=new_exp, env=env)
  rownames(resQTL) <- colnames(new_exp)
  colnames(resQTL) <- colnames(geno)
  write.table(resQTL, file=paste(location, gsub(".txt", "_QTL.txt", filename), sep=""))
  
  resInt <- apply(geno, 2, function(x, pheno, env){
    -log10(map.fast(x, pheno, env)$int)
  }, pheno=new_exp, env=env)
  rownames(resInt) <- colnames(new_exp)
  colnames(resInt) <- colnames(geno)
  write.table(resInt, file=paste(location, gsub(".txt", "_Int.txt", filename), sep=""))
  
  et <- proc.time()[3]
  cat(filename, " done in: ", et-st, "sec\n", sep="")
}
mapGenesByTU(filename = "genesByTU_chr1_norm_hf_cor.txt", geno, env)
