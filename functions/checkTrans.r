#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 29-01-2013
# first written: 24-01-2013
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
#show a strange trans on chr3!!!


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

checkTrans <- function(geno, env){
  location <- "C:/Arabidopsis Arrays/Data/OLD_unnormalized/new genes/"
  chrQTL <- vector("list", 6)
  chrMeansQTL <- NULL
  cnt <- 1
  for(filename in dir(location)[grepl("genes_by_chromosomes", dir(location))]){
    st <- proc.time()[3]
    new_exp <- t(read.table(paste("Data/OLD_unnormalized/new genes/", filename, sep=""), row.names=1, header=FALSE))
    chrQTL[[paste("chr", cnt, sep="")]] <- apply(geno, 2, function(x, pheno, env){
          -log10(map.fast(x, pheno, env)$qtl)
        }, pheno=new_exp, env=env)
    chrMeansQTL <- rbind(chrMeansQTL, colMeans(chrQTL[[paste("chr", cnt, sep="")]]))
    et <- proc.time()[3]
    cat("chr", cnt, " done in: ", et-st, "sec\n", sep="")
    cnt <- cnt + 1
  }
  chrQTL$chrMeansQTL <- chrMeansQTL
  return(chrQTL)
}

chrQTL <- checkTrans(geno, env)
save(chrQTL,  file="C:/Arabidopsis Arrays/Data/OLD_unnormalized/new genes/mapNewQTLgenes.Rdata")


load("C:/Arabidopsis Arrays/Data/OLD_unnormalized/new genes/mapNewQTLgenes.Rdata")
plotCheckTransP <- function(chrQTL){
  plot(c(1, 720), c(1, 70))
  for(chr in 1:5){
    for(x in 1:nrow(chrQTL[[paste("chr", chr, sep="")]])){
      points(1:716, chrQTL[[paste("chr", chr, sep="")]][x,], pch=20, col=chr)
    }
    axis(1, at=which.max(chrQTL$chrMeansQTL[chr,]), labels=which.max(chrQTL$chrMeansQTL[chr,]))
  }
}
plotCheckTransP(chrQTL)

plotCheckTransL <- function(chrQTL$chrMeansQTL){
  plot(c(1,720), c(1,max(chrQTL$chrMeansQTL)))
  for(chr in 1:5){
    points(chrQTL$chrMeansQTL[chr,], t='l', col=chr)
    axis(1, at=which.max(chrQTL$chrMeansQTL[chr,]), labels=which.max(chrQTL$chrMeansQTL[chr,]))
  }
}
plotCheckTransL(chrQTL)
