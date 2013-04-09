#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 16-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

getProbesOnChr <- function(ann_m, chr = 1){
  which(ann_m[,1] == chr)
}

getGeneStats <- function(filename, lodThreshold = 4, chrs = 1:5){
  probe_exp <- read.table(paste("Data/chr1/", filename, sep=""))
  gene_qtl <- read.table(paste("Data/chr1/", gsub(".txt", "_QTL.txt", filename), sep=""))
  sum_num <- NULL
  idmatrix <- vector("list", chrs)
  res <- NULL
  nprobes <- nrow(gene_qtl)
  for(chr in chrs){
    s <- 0
    ind <- NULL
    for(p in 1:nprobes){
      if(any(gene_qtl[p, getProbesOnChr(ann_m, chr)] >= lodThreshold)){
        s <- s+1
        ind <- c(ind, p)
      }
    }
    sum_num <- c(sum_num, s)
    if(is.null(ind)){
      idmatrix[[chr]] <- NA
    }else{
      idmatrix[[chr]] <- ind
    }
  }
  res$ratio <- sum_num / nprobes
  res$ind <- idmatrix
  res$means <- as.numeric(rowMeans(probe_exp[,17:ncol(probe_exp)]))
  res$nprobes <- nprobes
  return(res)
}



setwd("D:/Arabidopsis Arrays")
ann_m <- read.table("refined map/map.txt")
doc <- chr1
rawexp <- read.table(paste("Data/", doc, "/", fn_exp, sep=""), header=TRUE, row.names=1)
newexp <- rawexp[,17:164]
qtl <- read.table(paste("Data/", doc, "/", fn_qtl, sep=""))

getGeneStats <- function(newexp, qtl, lodThreshold = 4, chrs = 5){
  sum_num <- NULL
  idmatrix <- vector("list", chrs)
  res <- NULL
  nprobes <- nrow(qtl)
  for(chr in 1:chrs){
    s <- 0
    ind <- NULL
    for(p in 1:nprobes){
      if(any(qtl[p, getProbesOnChr(ann_m, chr)] >= lodThreshold)){
        s <- s+1
        ind <- c(ind, p)
      }
    }
    sum_num <- c(sum_num, s)
    if(is.null(ind)){
      idmatrix[[chr]] <- NA
    }else{
      idmatrix[[chr]] <- ind
    }
  }
  res$ratio <- sum_num / nprobes
  res$ind <- idmatrix
  res$means <- as.numeric(rowMeans(newexp[,17:ncol(newexp)]))
  res$nprobes <- nprobes
  return(res)
}
getGeneStats(newexp, qtl, lodThreshold = 4, chrs = 5)
