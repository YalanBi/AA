#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 29-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#******************************************************** this is the final version for reduced model mapping ^_^ *******************************************************#
setwd("D:/Arabidopsis Arrays")
geno <- read.table("refined map/genotypes.txt", sep="\t", row.names=1, header=TRUE)
menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]

map.fast <- function(geno, pheno, menvironment){
  res <- NULL
  modelinfo <- summary(aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno)))
  res$env   <- unlist(lapply(modelinfo, "[", 1, 5), use.names=TRUE)
  res$qtl   <- unlist(lapply(modelinfo, "[", 2, 5), use.names=TRUE)
  #res$int  <- unlist(lapply(modelinfo, "[", 3, 5), use.names=TRUE)
  res
}

mapGenotypes <- function(geno, menvironment, chr=1, verbose=TRUE){
  #env <- as.factor(menvironment)
  location <- paste0("Data/Raw/chr", chr, "_norm_hf_cor/")
  genenames <- gsub(".txt", "", dir(location)[grepl(".txt", dir(location)) & !grepl("QTL", dir(location)) & !grepl("SNP", dir(location))])
  
  for(filename in genenames){
    if(!file.exists(paste0(location, filename, "_QTL.txt"))){
      st <- proc.time()
      cat("Loading", filename, "\n")
      rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename), row.names=1, header=TRUE)
      if(nrow(rawexp) >= 4){
        newexp <- t(rawexp[ ,17:164])
     
        resQTL <- apply(geno, 2, function(x, pheno, env){
          -log10(map.fast(x, pheno, env)$qtl)
        }, pheno=newexp, env=menvironment)
        rownames(resQTL) <- rownames(rawexp)
        colnames(resQTL) <- colnames(geno)
        
        write.table(resQTL, file=paste0(location, filename, "_QTL.txt"))
        et <- proc.time()
        cat("Done with QTL mapping after:", (et-st)[3], "secs\n")
      }else{
        cat("Skipping", filename, "because it has to few probes\n")
      }
    }else{
      cat("Skipping", filename, "because it exists\n")
    }
  }
}

mapGenotypes(geno, menvironment, chr=1)
