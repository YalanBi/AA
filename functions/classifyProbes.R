#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 21-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

classifyProbes <- function(filename, expthreshold = 0.1){
  res <- getGeneStats(filename, lodThreshold = 4, chrs = 1:5) ############getGeneStats###########
  classInf <- NULL
  allprobes <- 1:res$nprobes
  qtlprobes <- unique(unlist(res$ind))
  expprobes <- which(res$means >= expmean)
  goodprobes <- unique(c(qtlprobes, expprobes))
  badprobes <- which(!allprobes %in% goodprobes)
  classInf$badP <- badprobes
  classInf$goodP <- goodprobes
  classInf$qtlP <- qtlprobes
  classInf$expP <- expprobes
  return(classInf)
}




#New version
setwd("C:/Arabidopsis Arrays")
ann_m <- read.table("refined map/map.txt")

getProbesOnChr <- function(ann_m, chr = 1){
  which(ann_m[,1] == chr)
}

getGeneStats <- function(filename, lodThreshold = 4, chrs = 5){
  qtl <- read.table(paste("Data/chr1/", gsub(".txt", "_QTL.txt", filename), sep=""))
  nprobes <- nrow(qtl)
  sum_num <- NULL
  idmatrix <- vector("list", chrs)
  res <- NULL
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
  res$nprobes <- nprobes
  return(res)
}

expprobes <- function(newexp, ratio = 0.1){
  expprb_ind <- NULL
  for(a in 1:nrow(newexp)){
    if(sum(newexp[a,] >= 4.5)/ncol(newexp) >= ratio){
      expprb_ind <- c(expprb_ind, a)
    }
  }
  return(expprb_ind)
}

classifyProbes <- function(filename, ratio){
  rawexp <- read.table(paste("Data/", doc, "/", filename, sep=""), header=TRUE, row.names=1)
  newexp <- rawexp[,17:164]
  res <- getGeneStats(filename, lodThreshold = 4, chrs = 5)
  classInf <- NULL
  allprobes <- 1:res$nprobes
  qtlprobes <- unique(unlist(res$ind))
  expprobes <- expprobes(newexp, ratio)
  goodprobes <- unique(c(qtlprobes, expprobes))
  badprobes <- which(!allprobes %in% goodprobes)
  classInf$badP <- badprobes
  classInf$goodP <- goodprobes
  classInf$qtlP <- qtlprobes
  classInf$expP <- expprobes
  return(classInf)
}

classifyProbes(filename, ratio = 0.5)
