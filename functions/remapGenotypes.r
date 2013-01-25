#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 24-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("C:/Arabidopsis Arrays")
geno_new <- read.table("refined map/genotypes.txt", row.names=1)
geno <- read.table("Data/Old/genotypes_n.txt", row.names=1)
ann <- read.table("Data/ann_env.txt", colClasses="character")[,2]
new_exp <- t(read.table("genes_by_chromosomes.txt",row.names=1,header=FALSE))


map.fast <- function(x, geno, pheno, menvironment){
	res <- NULL
  envv <- as.factor(menvironment)
	models  <- aov(as.matrix(pheno) ~ envv + as.numeric(geno[,x]) +  envv:as.numeric(geno[,x]))
	modelinfo <- summary(models)
	res$env <- unlist(lapply(modelinfo,"[",1,5),use.names=T)
	res$qtl <- unlist(lapply(modelinfo,"[",2,5),use.names=T)
  res$int <- unlist(lapply(modelinfo,"[",3,5),use.names=T)
	res
}

qtlmatrix <- NULL
intmatrix <- NULL
geno <- geno[rownames(geno_new),]
for(m in 35){
  res <- map.fast(m, geno,new_exp, menvironment)
  qtlmatrix <- cbind(qtlmatrix,res$qtl)
  intmatrix <- cbind(intmatrix,res$int)
  cat(m,"\n")
}
