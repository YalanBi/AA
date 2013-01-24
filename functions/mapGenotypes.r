#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 16-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("C:/Arabidopsis Arrays")
geno <- read.table("refined map/genotypes.txt",sep="\t",row.names=1,header=TRUE)
menvironment <- read.table("Data/ann_env.txt",sep="\t")[,2]

map.fast <- function(geno, pheno, menvironment){
  res <- NULL
  modelinfo <- summary(aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno)))
  res$env   <- unlist(lapply(modelinfo,"[",1,5),use.names=T)
  res$qtl   <- unlist(lapply(modelinfo,"[",2,5),use.names=T)
  #res$int   <- unlist(lapply(modelinfo,"[",3,5),use.names=T)
  res
}

mapGenotypes <- function(geno, menvironment, chr=1){
  env <- as.factor(menvironment)
  for(filename in dir(paste("Data/chr", chr, sep=""))[grepl(".txt",dir(paste("Data/chr", chr, sep="")))]){
    if(!file.exists(paste("Data/chr", chr, "/", gsub(".txt","_QTL.txt",filename), sep="")) && !grepl("_QTL",filename)){
      st <- proc.time()
      cat("Loading", filename,"\n")
      Pheno <- read.table(paste("Data/chr", chr, "/", filename, sep=""),row.names=1, header=TRUE)
      if(nrow(Pheno) < 4){
        cat("Skipping", filename," because it has to few probes\n")
      }else{
        pheno <- t(Pheno[,16:ncol(Pheno)])#####col numbers are changed!!!#####
     
        resQTL <- apply(geno, 2, function(x, pheno, env){
          -log10(map.fast(x, pheno, env)$qtl)
        }, pheno=pheno, env=env)
        rownames(resQTL) <- colnames(pheno)
        colnames(resQTL) = colnames(geno)
        write.table(resQTL, file=paste("Data/chr", chr, "/", gsub(".txt","_QTL.txt",filename), sep=""))
        et <- proc.time()
        cat("Done with QTL mapping after:",(et-st)[3],"secs\n")
      }
    }else{
      cat("Skipping", filename," because it exists\n")
    }
  }
}

mapGenotypes(geno, menvironment, chr=1)
