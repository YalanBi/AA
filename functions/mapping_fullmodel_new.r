#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 17-06-2013
# first written: 22-03-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("D:/Arabidopsis Arrays")
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)
menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]


map.fast <- function(geno, pheno, menvironment){
  res <- NULL
  models    <- aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno) + as.factor(menvironment):as.numeric(geno))
  modelinfo <- summary(models)
  res$env   <- unlist(lapply(modelinfo, "[", 1, 5), use.names = T)
  res$qtl   <- unlist(lapply(modelinfo, "[", 2, 5), use.names = T)
  res$int   <- unlist(lapply(modelinfo, "[", 3, 5), use.names = T)
  res$eff   <- unlist(models$coefficients[5, ])
  res
}


fullModelMapping <- function(filename, geno, menvironment, P = -1, verbose = FALSE){
  st <- proc.time()

  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  
  #if P is positive, check nprobes
  if(nrow(rawexp) >= P){
    if(verbose) cat(filename, "is >=", P, "probes, start mapping\n")
    pheno <- t(rawexp[,17:164])

    resList <- apply(geno, 2, function(x, pheno, env){
      result <- map.fast(x, pheno, menvironment)
      return(list(-log10(result$env), sign(result$eff) * -log10(result$qtl), -log10(result$int)))
    }, pheno = pheno, env = menvironment)

    qtls <- lapply(resList, "[[", 2)
    resQTL <- matrix(unlist(qtls), length(qtls[[1]]), length(qtls))
    rownames(resQTL) <- rownames(rawexp)
    colnames(resQTL) <- colnames(geno)
    write.table(resQTL, file = paste0("Data/FullModel/chr", chr, "_norm_hf_cor_FM/", filename, "_FM_QTL.txt"))

    envs <- lapply(resList, "[[", 1)
    resEnv <- matrix(unlist(envs), length(envs[[1]]), length(envs))
    rownames(resEnv) <- rownames(rawexp)
    colnames(resEnv) <- colnames(geno)
    write.table(resEnv, file = paste0("Data/FullModel/chr", chr, "_norm_hf_cor_FM/", filename, "_FM_Env.txt"))

    ints <- lapply(resList, "[[", 3)
    resInt <- matrix(unlist(ints), length(ints[[1]]), length(ints))
    rownames(resInt) <- rownames(rawexp)
    colnames(resInt) <- colnames(geno)
    write.table(resInt, file = paste0("Data/FullModel/chr", chr, "_norm_hf_cor_FM/", filename, "_FM_Int.txt"))

    et <- proc.time()
    if(verbose) cat(filename, "done after:", (et-st)[3], "secs\n")
  } else if(verbose) cat("skip", filename, ", because less than", P, "probes\n")
}


for(chr in 1:5){
  st <- proc.time()
  genenames <- gsub(".txt", "", dir(paste0("Data/Raw/chr", chr, "_norm_hf_cor/"))[!grepl("SNP", dir(paste0("Data/Raw/chr", chr, "_norm_hf_cor/")))])
  #filename = "AT1G01010"
  for(filename in genenames){
    if(!file.exists(paste0("Data/FullModel/chr", chr, "_norm_hf_cor_FM/", filename, "_FM_QTL.txt"))){
      fullModelMapping(filename, geno, menvironment, P = 4)
    } else cat("Skipping", filename, ", because it exists\n")
  }
  et <- proc.time()
  cat("chr", chr, "is done with QTL mapping, after:", (et-st)[3], "secs\n")
}

