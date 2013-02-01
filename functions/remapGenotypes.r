#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 29-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#map by each chr
setwd("C:/Arabidopsis Arrays")
geno <- read.table("refined map/genotypes.txt", row.names=1)
#geno <- read.table("Data/Old/genotypes_n.txt", row.names=1)
menvironment <- read.table("Data/ann_env.txt",sep="\t")[,2]

new_exp <- t(read.table("Data/OLD_unnormalized/new genes/genes_by_chromosomes2.txt",row.names=1,header=FALSE))
#new_exp <- t(read.table("Data/OLD_unnormalized/Chr2/At2G01190.txt",row.names=1,header=FALSE)[,16:163])
#new_exp <- new_exp[,1:200]

map.fast <- function(x, geno, pheno, menvironment){
  res <- NULL
  #envv <- as.factor(menvironment)
  models  <- aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno[,x]) +  as.factor(menvironment):as.numeric(geno[,x]))
  modelinfo <- summary(models)
  res$env <- unlist(lapply(modelinfo,"[",1,5),use.names=T)
  res$qtl <- unlist(lapply(modelinfo,"[",2,5),use.names=T)
  res$int <- unlist(lapply(modelinfo,"[",3,5),use.names=T)
  res
}

qtlmatrix <- NULL
intmatrix <- NULL
for(m in 1:ncol(geno)){
  res <- map.fast(m, geno,new_exp, menvironment)
  qtlmatrix <- cbind(qtlmatrix,res$qtl)
  intmatrix <- cbind(intmatrix,res$int)
  cat(m,"\n")
}
image(t(-log10(qtlmatrix))>5)
#show a strange trans on chr3 before normalization!!!


#unbalanced individuals of env 2 on the 326th marker cause this trans!
hist(as.numeric(menvironment)+(geno[,326]-1.5)/10, breaks=20)



#map.fast for 5 chrs
setwd("C:/Arabidopsis Arrays")
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

mapNewQTLGenes <- function(filename = "genes_by_chromosomes1_normalized.txt", geno, env){
  location <- "C:/Arabidopsis Arrays/Data/"
  st <- proc.time()[3]
  new_exp <- t(read.table(paste("Data/", filename, sep=""), row.names=1, header=FALSE))
  res <- -log10(apply(geno, 2, function(x, pheno, env){
    map.fast(x, pheno, env)$qtl
  }, pheno=new_exp, env=env))
  write.table(res, file=paste(location, gsub(".txt", "_QTL.txt", filename), sep=""))
  et <- proc.time()[3]
  cat(filename, " done in: ", et-st, "sec\n", sep="")
}
#mapNewQTLGenes(filename = "genes_by_chromosomes1_normalized.txt", geno, env)


image(t(res)>4)



mapNewGOODGenes <- function(filename = "genes_summarized_1_normalized.txt", geno, env){
  location <- "C:/Arabidopsis Arrays/Data/"
  st <- proc.time()[3]
  new_exp <- t(read.table(paste("Data/", filename, sep=""), row.names=1, header=FALSE))
  res <- NULL
  res <- apply(geno, 2, function(x, pheno, env){
    -log10(map.fast(x, pheno, env)$qtl)
  }, pheno=new_exp, env=env)
  write.table(res, file=paste(location, gsub(".txt", "_QTL.txt", filename), sep=""))
  et <- proc.time()[3]
  cat(filename, " done in: ", et-st, "sec\n", sep="")
}
#mapNewGOODGenes(filename = "genes_summarized_1_normalized.txt", geno, env)





#version 0_0 for unnormalized data
#map.fast for 5 chrs
setwd("C:/Arabidopsis Arrays")
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

mapNewQTLGenes <- function(geno, env){
  location <- "C:/Arabidopsis Arrays/Data/OLD_unnormalized/new genes/"
  cnt <- 1
  for(filename in dir(location)[grepl("genes_by_chromosomes", dir(location))]){
    st <- proc.time()[3]
    new_exp <- t(read.table(paste("Data/OLD_unnormalized/new genes/", filename, sep=""), row.names=1, header=FALSE))
    res <- NULL
    res <- apply(geno, 2, function(x, pheno, env){
      -log10(map.fast(x, pheno, env)$qtl)
    }, pheno=new_exp, env=env)
    write.table(res, file=paste("Data/OLD_unnormalized/new genes/genes_by_chromosomes", cnt, "_QTL.txt", sep=""))
    et <- proc.time()[3]
    cat("chr", cnt, " done in: ", et-st, "sec\n", sep="")
    cnt <- cnt + 1
  }
}
mapNewQTLGenes(geno, env)


image(t(res)>4)



mapNewGOODGenes <- function(geno, env){
  location <- "C:/Arabidopsis Arrays/Data/OLD_unnormalized/new genes/"
  cnt <- 1
  for(filename in dir(location)[grepl("genes_summarized", dir(location))]){
    st <- proc.time()[3]
    new_exp <- t(read.table(paste("Data/OLD_unnormalized/new genes/", filename, sep=""), row.names=1, header=FALSE))
    res <- NULL
    res <- apply(geno, 2, function(x, pheno, env){
      -log10(map.fast(x, pheno, env)$qtl)
    }, pheno=new_exp, env=env)
    write.table(res, file=paste("Data/OLD_unnormalized/new genes/genes_summarized_chr", cnt, "_QTL.txt", sep=""))
    et <- proc.time()[3]
    cat("chr", cnt, " done in: ", et-st, "sec\n", sep="")
    cnt <- cnt + 1
  }
}

mapNewGOODGenes(geno, env)
