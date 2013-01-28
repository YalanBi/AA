

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

checkTrans <- function(geno, env){
  location <- "C:/Arabidopsis Arrays/Data/new genes/"
  chrMeansQTL <- NULL
  cnt <- 1
  plot(c(1,750), c(1,8))
  for(filename in dir(location)[grepl("genes_by_chromosomes", dir(location))]){
    st <- proc.time()[3]
    new_exp <- t(read.table(paste("Data/new genes/", filename, sep=""), row.names=1, header=FALSE))
    chrMeansQTL <- colMeans(apply(geno, 2, function(x, pheno, env){
          -log10(map.fast(x, pheno, env)$qtl)
        }, pheno=new_exp, env=env))
      points(chrMeansQTL, t='l', col=cnt)
    et <- proc.time()[3]
    cat("chr", cnt, " done in: ", et-st, "sec\n", sep="")
    cnt <- cnt + 1
  }
}

checkTrans(geno, env)