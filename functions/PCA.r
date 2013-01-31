#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 29-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("C:/Arabidopsis Arrays")
geno <- read.table("refined map/genotypes.txt",sep="\t",row.names=1,header=TRUE)
menvironment <- read.table("Data/ann_env.txt",sep="\t")[,2]

newP <- read.table("Data/new10kProbes.txt", row.names=1, header=TRUE)
pcares <- prcomp(newP[,17:ncol(newP)])
write.table(pcares[[2]], file="Data/pcaRes.txt")
plot(pcares[[2]][,1],t='h')

#map only the 1st PC
map.fast <- function(x, geno, pheno, menvironment){
  res <- NULL
  modelinfo <- summary(aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno[,x]) +  as.factor(menvironment):as.numeric(geno[,x])))
  res$env <- unlist(lapply(modelinfo,"[",1,5),use.names=T)
  res$qtl <- unlist(lapply(modelinfo,"[",2,5),use.names=T)
  res$int <- unlist(lapply(modelinfo,"[",3,5),use.names=T)
  return(res)
}

pcaqtl <- NULL
for(m in 1:ncol(geno)){
  pcaqtl <- c(pcaqtl, -log10(map.fast(m, geno, pcares[[2]][,1], menvironment)$qtl))
  cat(m,"\n")
}

which.max(pcaqtl)
#[1] 59
which(pcaqtl>3)
#[1]  47  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67  71  72  73  74  75  76  77 321 322 323 324 325 326
hist(pcaqtl,breaks=25)
plot(pcaqtl,t='l')



#map PC1-148
map.fast <- function(geno, pheno, menvironment){
  res <- NULL
  modelinfo <- summary(aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno) + as.factor(menvironment):as.numeric(geno)))
  res$env   <- unlist(lapply(modelinfo,"[",1,5), use.names=T)
  res$qtl   <- unlist(lapply(modelinfo,"[",2,5), use.names=T)
  res$int   <- unlist(lapply(modelinfo,"[",3,5), use.names=T)
  return(res)
}

pcaqtlm <- res <- apply(geno, 2, function(x, pheno, env){
      -log10(map.fast(x, pheno, env)$qtl)
    }, pheno=pcares[[2]], env=menvironment)
rownames(pcaqtlm) <- gsub(" Response ", "qtl_", rownames(pcaqtlm))
write.table(pcaqtlm, file="Data/pcaQTL.txt")
