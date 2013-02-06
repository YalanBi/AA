#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 29-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("C:/Arabidopsis Arrays")
ann_m <- read.table("refined map/map.txt")

getProbesOnChr <- function(ann_m, chr = 1){
  which(ann_m[,1] == chr)
}


createGenes <- function(chromosome = 1, empty = TRUE){  #TODO: Per chromosome we need a new file !!!
  location <- paste("C:/Arabidopsis Arrays/Data/chr",chromosome,"_norm_hf_cor/",sep="")
  if(empty) cat("", file=paste("Data/genes_by_chromosomes", chromosome, "_norm_hf_cor.txt", sep=""))
  load(paste(location,"Classification_chr",chromosome,"_norm_hf_cor.Rdata",sep=""))
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
      probesa4 <- good[which(apply(qtl_data[good,marker_ids],1,max) >= 4)] #lodThreshold=4
      if(length(probesa4) >= 4){
        #Create a new gene !
        newgenename <- gsub(".txt", paste("_",chr,sep=""), fname)
        expValues <- colMeans(exp_data[probesa4, 17:ncol(exp_data)])
        cat(newgenename, expValues, "\n", file=paste("Data/genes_by_chromosomes", chromosome, "_norm_hf_cor.txt", sep=""), append=TRUE)
        cnt <- cnt+1
      }
    }
    cat("Finished with", fname, "created:", cnt, "Genes\n")
  }
}
createGenes(chromosome = 2)


createSingleGenes <- function(chromosome = 1, empty = TRUE){
  location <- paste("C:/Arabidopsis Arrays/Data/chr", chromosome, "_norm_hf_cor/", sep="")
  if(empty) cat("", file=paste("Data/genes_summarized_", chromosome, "_norm_hf_cor.txt", sep=""))
  load(paste(location,"Classification_chr",chromosome,"_norm_hf_cor.Rdata",sep=""))
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
    if(length(good) > 0){
      expValues <- colMeans(exp_data[good, 17:ncol(exp_data)])
      newgenename <- gsub(".txt", "", fname)
      cat(newgenename, expValues, "\n", file=paste("Data/genes_summarized_", chromosome, "_norm_hf_cor.txt", sep=""), append=TRUE)
      cat("Finished with", fname, "\n")
    }else{
      cat(fname, "no good probes\n")
    }
  }
}
createSingleGenes(chromosome = 2)


for(n in 1:5){
  if(file.exists(paste("Data/genes_by_chromosomes", n, "_norm_hf_cor.txt", sep=""))){
    cat(paste("Data/genes_by_chromosomes", n, "_norm_hf_cor.txt", sep=""), "exist\n")
  }
  else{
    createGenes(chromosome = n)
  }
  if(file.exists(paste("Data/genes_summarized_", n, "_norm_hf_cor.txt", sep=""))){
    cat(paste("Data/genes_summarized_", n, "_norm_hf_cor.txt", sep=""), "exist\n")
  }
  else{
    createSingleGenes(chromosome = n)
  }
}




#0_0 version for unnormalized data
setwd("C:/Arabidopsis Arrays")
ann_m <- read.table("refined map/map.txt")

getProbesOnChr <- function(ann_m, chr = 1){
  which(ann_m[,1] == chr)
}


createGenes <- function(chromosome = 1, empty = TRUE){  #TODO: Per chromosome we need a new file !!!
  location <- paste("C:/Arabidopsis Arrays/Data/OLD_unnormalized/chr", chromosome, "/", sep="")
  if(empty) cat("", file=paste("C:/Arabidopsis Arrays/Data/OLD_unnormalized/new genes/genes_by_chromosomes", chromosome, ".txt", sep=""))
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
      probesa4 <- good[which(apply(qtl_data[good,marker_ids],1,max) >= 4)] #lodThreshold=4
      if(length(probesa4) >= 4){
        #Create a new gene !
        newgenename <- gsub(".txt", paste("_",chr,sep=""), fname)
        expValues <- colMeans(exp_data[probesa4, 16:ncol(exp_data)])
        cat(newgenename, expValues, "\n", file=paste("C:/Arabidopsis Arrays/Data/OLD_unnormalized/new genes/genes_by_chromosomes", chromosome, ".txt", sep=""), append=TRUE)
        cnt <- cnt+1
      }
    }
    cat("Finished with", fname, "created:", cnt, "Genes\n")
  }
}
createGenes(chromosome = 1)


createSingleGenes <- function(chromosome = 1, empty = TRUE){
  location <- paste("C:/Arabidopsis Arrays/Data/OLD_unnormalized/chr", chromosome, "/", sep="")
  if(empty) cat("", file=paste("C:/Arabidopsis Arrays/Data/OLD_unnormalized/new genes/genes_summarized_chr", chromosome, ".txt", sep=""))
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
    if(length(good) > 0){
      expValues <- colMeans(exp_data[good, 16:ncol(exp_data)])
      newgenename <- gsub(".txt", "", fname)
      cat(newgenename, expValues, "\n", file=paste("C:/Arabidopsis Arrays/Data/OLD_unnormalized/new genes/genes_summarized_chr", chromosome, ".txt", sep=""), append=TRUE)
      cat("Finished with", fname, "\n")
    }else{
      cat(fname, "no good probes\n")
    }
  }
}
createSingleGenes(chromosome = 1)
