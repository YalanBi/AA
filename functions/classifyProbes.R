#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 29-01-2013
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
setwd("D:/Arabidopsis Arrays")
ann_m <- read.table("refined map/map.txt")

getProbesOnChr <- function(ann_m, chr = 1){
  which(ann_m[,1] == chr)
}

getGeneStats <- function(filename, lodThreshold = 4, chrs = 5){
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
    #Here we know which probes we're going to summarize
    
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

#Simplify, without idmatrix
getGeneStatsS <- function(filename, lodThreshold = 4, chrs = 5){
  #st <- proc.time()[3]
  qtl <- read.table(filename)
  nprobes <- nrow(qtl)
  res <- NULL
  ind <- NULL
  for(p in 1:nprobes){
    if(any(qtl[p,] >= lodThreshold)){
      ind <- c(ind, p)
    }
  }
  if(!is.null(ind)){
    res$qtlind <- ind
  }
  else{
    res$qtlind <- NULL
  }
   res$nprobes <- nprobes
   #et <- proc.time()[3]
   #cat(et-st, "\n")
   return(res)
}
#getGeneStatsS("C:\\Arabidopsis Arrays\\Data\\chr2_normalized\\AT2G01060_QTL.txt", lodThreshold = 4, chrs = 5)

#2 variables to control
expprobes_ratio <- function(newexp, expThreshold=4.5, ratio = 0.1){
  expprb_ind <- NULL
  for(a in 1:nrow(newexp)){
    if(sum(newexp[a,] >= expThreshold)/ncol(newexp) >= ratio){
      expprb_ind <- c(expprb_ind, a)
    }
  }
  return(expprb_ind)
}

expprobes_mean <- function(newexp, expThreshold=5){
  expprb_ind <- NULL
  for(a in 1:nrow(newexp)){
    if(mean(as.numeric(newexp[a,])) >= expThreshold){
      expprb_ind <- c(expprb_ind, a)
    }
  }
  return(expprb_ind)
}

classifyProbes <- function(filename, lodThreshold = 4, expThreshold = 5, chrs = 5){
  rawexp <- read.table(file=gsub("_QTL", "", filename), header=TRUE, row.names=1)
  newexp <- rawexp[,17:164]
  #res <- getGeneStats(filename, lodThreshold = lodThreshold, chrs = chrs)#with idmatrix
  res <- getGeneStatsS(filename, lodThreshold = lodThreshold, chrs = chrs)
  classInf <- NULL
  allprobes <- 1:res$nprobes
  #qtlprobes <- unique(unlist(res$ind))
  qtlprobes <- res$qtlind
  expprobes <- expprobes_mean(newexp, expThreshold=expThreshold)
  intronprobes <- which(! allprobes %in% grep("tu", rawexp[,"tu"]))
  goodprobes <- unique(c(qtlprobes, expprobes, intronprobes))# A intron probe is always GOOD !
  classInf$badP <- which(!allprobes %in% goodprobes)
  classInf$goodP <- goodprobes[order(goodprobes, decreasing = FALSE)]
  #classInf$qtlP <- qtlprobes[order(qtlprobes, decreasing = FALSE)]
  classInf$qtlP <- qtlprobes
  classInf$expP <- expprobes
  classInf$introP <- intronprobes
  return(classInf)
}

setwd("C:\\Arabidopsis Arrays\\Data\\chr1_norm_hf_cor")
res <- list()
for(x in dir()[grepl("_QTL",dir())]){
  st <- proc.time()
  res[[x]] <- classifyProbes(filename = x) #Please note: Filename = QTL file !!
  et <- proc.time()
  cat(x, "done after:", (et-st)[3], "secs\n")
}
save(res,  file="Classification_chr1_norm_hf_cor.Rdata")
