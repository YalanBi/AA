#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 24-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("C:/Arabidopsis Arrays")
ann_m <- read.table("refined map/map.txt")

getProbesOnChr <- function(ann_m, chr = 1){
  which(ann_m[,1] == chr)
}


createGenes <- function(chromosome = 1, empty = TRUE){  #TODO: Per chromosome we need a new file !!!
  location <- paste("C:/Arabidopsis Arrays/Data/chr",chromosome,"/",sep="")
  if(empty) cat("", file="genes_by_chromosomes.txt")
  load(paste(location,"Classification_chr",chromosome,".Rdata",sep=""))
  ann_m <- read.table("C:/Arabidopsis Arrays/refined map/map.txt")

  for(filename in dir(location)[grepl("_QTL",dir(location))]){
    qtl_data <- read.table(paste(location,filename,sep=""))
    fname <- gsub("_QTL.txt",".txt", filename)
    exp_data <- read.table(paste(location,fname,sep=""),header=T,row.names=1)
    class <- res[[filename]]
    if(length(which(class$goodP %in% class$introP)) > 0){
      good <-  class$goodP[- which(class$goodP %in% class$introP)]
    }else{
      good <-  class$goodP
    }
    cnt <- 0
    for(chr in 1:5){
      marker_ids <- getProbesOnChr(ann_m, chr)
      probesa5 <- good[which(apply(qtl_data[good,marker_ids],1,max) > 5)]
      if(length(probesa5) > 4){
        #Create a new gene !
        newgenename <- gsub(".txt", paste("_",chr,sep=""), fname)
        expValues <- colMeans(exp_data[probesa5, 16:ncol(exp_data)])
        cat(newgenename, expValues, "\n", file="genes_by_chromosomes.txt", append=TRUE)
        cnt <- cnt+1
      }
    }
    cat("Finished with", fname, "created:", cnt, "Genes\n")
  }
}


createSingleGenes <- function(chromosome = 1, empty = TRUE){
  location <- paste("C:/Arabidopsis Arrays/Data/chr",chromosome,"/",sep="")
  if(empty) cat("", file=paste("genes_summarized_",chromosome,".txt",sep=""))
  load(paste(location,"Classification_chr",chromosome,".Rdata",sep=""))
  ann_m <- read.table("C:/Arabidopsis Arrays/refined map/map.txt")

  for(filename in dir(location)[grepl("_QTL",dir(location))]){
    qtl_data <- read.table(paste(location,filename,sep=""))
    fname <- gsub("_QTL.txt",".txt", filename)
    exp_data <- read.table(paste(location,fname,sep=""),header=T,row.names=1)
    class <- res[[filename]]
    if(length(which(class$goodP %in% class$introP)) > 0){
      good <-  class$goodP[- which(class$goodP %in% class$introP)]
    }else{
      good <-  class$goodP
    }
    expValues <- colMeans(exp_data[good, 16:ncol(exp_data)])
    newgenename <- gsub(".txt", "", fname)
    cat(newgenename, expValues, "\n", file="genes_summarized.txt", append=TRUE)
    cat("Finished with", fname, "\n")
  }
}




setwd("C:/Arabidopsis Arrays")
#geno <- read.table("refined map/genotypes.txt", row.names=1)
geno <- read.table("Data/Old/genotypes_n.txt", row.names=1)
ann <- read.table("Data/ann_env.txt", colClasses="character")
menvironment <- ann[,2]
new_exp <- t(read.table("genes_by_chromosomes.txt",row.names=1,header=FALSE))

new_exp <- new_exp[,1:400]

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
for(m in 1:ncol(geno)){
  res <- map.fast(m, geno,new_exp, menvironment)
  qtlmatrix <- cbind(qtlmatrix,res$qtl)
  intmatrix <- cbind(intmatrix,res$int)
  cat(m,"\n")
}
