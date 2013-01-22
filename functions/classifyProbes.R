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

getGeneStats <- function(filename, lodThreshold = 5, chrs = 5){
  qtl <- read.table(filename)
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
    if(!is.null(ind)){
      idmatrix[[chr]] <- ind
    }
  }
  res$ratio <- sum_num / nprobes
  res$ind <- idmatrix
  res$nprobes <- nprobes
  return(res)
}

expprobes <- function(newexp, expThreshold=4.5, ratio = 0.1){
  expprb_ind <- NULL
  for(a in 1:nrow(newexp)){
    if(sum(newexp[a,] >= expThreshold)/ncol(newexp) >= ratio){
      expprb_ind <- c(expprb_ind, a)
    }
  }
  return(expprb_ind)
}

expprobes_mean <- function(newexp, expThreshold=4.5){
  expprb_ind <- NULL
  for(a in 1:nrow(newexp)){
    if(mean(as.numeric(newexp[a,])) >= expThreshold){
      expprb_ind <- c(expprb_ind, a)
    }
  }
  return(expprb_ind)
}

classifyProbes <- function(filename, lodThreshold = 5, expThreshold = 5, ratio = 0.1, chrs = 5){
  rawexp <- read.table(file=gsub("_QTL", "", filename), header=TRUE, row.names=1)
  newexp <- rawexp[,17:164]
  res <- getGeneStats(filename, lodThreshold = lodThreshold, chrs = chrs)
  classInf <- NULL
  allprobes <- 1:res$nprobes
  qtlprobes <- unique(unlist(res$ind))
  expprobes <- expprobes_mean(newexp, expThreshold=expThreshold)
  intronprobes <- which(! allprobes %in% grep("tu", rawexp[,"tu"]))
  goodprobes <- unique(c(qtlprobes, expprobes, intronprobes))# A intron probe is always GOOD !
  classInf$badP <- which(!allprobes %in% goodprobes)
  classInf$goodP <- goodprobes[order(goodprobes, decreasing = FALSE)]
  classInf$qtlP <- qtlprobes
  classInf$expP <- expprobes
  classInf$introP <- intronprobes
  return(classInf)
}

setwd("C:\\Arabidopsis Arrays\\Data\\chr4")
res <- list()
for(x in dir()[grepl("_QTL",dir())]){
  st <- proc.time()
  res[[x]] <- classifyProbes(filename = x) #Please note: Filename = QTL file !!
  et <- proc.time()
  cat(x, "done after:", (et-st)[3], "secs\n")
}
save(res,  file="Classification_chr4.Rdata")




#only bad probes
classifyBadP <- function(filename, lodThreshold = 5, expThreshold = 5, ratio = 0.1, chrs = 5){
  rawexp <- read.table(file=gsub("_QTL", "", filename), header=TRUE, row.names=1)
  newexp <- rawexp[,17:164]
  res <- getGeneStats(filename, lodThreshold = lodThreshold, chrs = chrs)
  classInf <- NULL
  allprobes <- 1:res$nprobes
  qtlprobes <- unique(unlist(res$ind))
  expprobes <- expprobes_mean(newexp, expThreshold=expThreshold)
  intronprobes <- which(! allprobes %in% grep("tu", rawexp[,"tu"]))
  goodprobes <- unique(c(qtlprobes, expprobes, intronprobes))# A intron probe is always GOOD !
  classInf <- which(!allprobes %in% goodprobes)
  return(classInf)
}
