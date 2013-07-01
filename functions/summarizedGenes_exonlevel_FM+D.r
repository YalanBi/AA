#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 01-07-2013
# first written: 19-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#****************************************** this is the final version for full model mapping with directions at exon level! ^_^ ******************************************#
#************************************** one probe for one exon, and the order of RILs is as same as their initial order (no change) **************************************#
#********************************************* this is used for PERMUTAION, CIS-TRANS PLOT and MAIN/CONSISTENT eQTL and Int! *********************************************#

setwd("D:/Arabidopsis Arrays")
geno <- read.table("refined map/genotypes.txt", sep="\t", row.names=1, header=TRUE)
menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]

#function for full model mapping
map.fast <- function(geno, pheno, menvironment){
  res <- NULL
  models    <- aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno) + as.factor(menvironment):as.numeric(geno))
  modelinfo <- summary(models)
  res$env   <- unlist(lapply(modelinfo, "[", 1, 5), use.names=TRUE)
  res$qtl   <- unlist(lapply(modelinfo, "[", 2, 5), use.names=TRUE)
  res$int   <- unlist(lapply(modelinfo, "[", 3, 5), use.names=TRUE)
  res$eff   <- unlist(models$coefficients[5, ])
  res
}

#full model mapping with summarized expressed genes
fullModelMapping_exonLevel <- function(chr=1, geno, menvironment){
  st <- proc.time()
  cat("Loading chr", chr, "...\n")
  #the 1st col is bp!!!
  rawexp <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_1tu1probe.txt"), row.names=1, header=TRUE)[ ,2:149]
  pheno <- t(rawexp)
  
  resEnv <- NULL
  resQTL <- NULL
  resInt <- NULL
  
  for(m in 1:ncol(geno)){
    resList <- map.fast(geno[ ,m], pheno=pheno, menvironment)
    resEnv <- cbind(resEnv, -log10(resList$env))
    resQTL <- cbind(resQTL, -log10(resList$qtl) * sign(resList$eff))
    resInt <- cbind(resInt, -log10(resList$int))
  }
  
  rownames(resEnv) <- rownames(rawexp)
  colnames(resEnv) <- colnames(geno)
  write.table(resEnv, file=paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_Env.txt"), sep="\t")
  
  rownames(resQTL) <- rownames(rawexp)
  colnames(resQTL) <- colnames(geno)
  write.table(resQTL, file=paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_QTL.txt"), sep="\t")
  
  rownames(resInt) <- rownames(rawexp)
  colnames(resInt) <- colnames(geno)
  write.table(resInt, file=paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_Int.txt"), sep="\t")
  
  et <- proc.time()
  cat("Done with full model mapping after:", (et-st)[3], "secs\n")
}

#full model mapping, 5 chr
for(chr in 1:5){
  fullModelMapping_exonLevel(chr = chr, geno, menvironment)
}


#another function for full model mapping, same results
fullModelMapping_exonLevel <- function(chr, geno, menvironment){
  st <- proc.time()
  cat("Loading chr", chr, "...\n")
  rawexp <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_1tu1probe.txt"), row.names=1, header=TRUE)[ ,2:149]
  pheno <- t(rawexp)
  
  resList <- apply(geno, 2, function(x, pheno, env){
    result <- map.fast(x, pheno, menvironment)
    return(list(-log10(result$env), sign(result$eff) * -log10(result$qtl), -log10(result$int)))
  }, pheno=pheno, env=menvironment)
  
  qtls <- lapply(resList, "[[", 2)
  resQTL <- matrix(unlist(qtls), length(qtls[[1]]), length(qtls))#nrow(resQTL)=length(qtls[[1]]); ncol(resQTL)=length(qtls)
  rownames(resQTL) <- rownames(rawexp)
  colnames(resQTL) <- colnames(geno)
  write.table(resQTL, file = paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_QTL.txt"))
  
  envs <- lapply(resList, "[[", 1)
  resEnv <- matrix(unlist(envs), length(envs[[1]]), length(envs))
  rownames(resEnv) <- rownames(rawexp)
  colnames(resEnv) <- colnames(geno)
  write.table(resEnv, file = paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_Env.txt"))
  
  ints <- lapply(resList, "[[", 3)
  resInt <- matrix(unlist(ints), length(ints[[1]]), length(ints))
  rownames(resInt) <- rownames(rawexp)
  colnames(resInt) <- colnames(geno)
  write.table(resInt, file = paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_Int.txt"))
  
  et <- proc.time()
  cat("Done with full model mapping after:", (et-st)[3], "secs\n")
}
