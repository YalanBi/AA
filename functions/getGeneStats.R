#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 16-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

getProbesOnChr <- function(ann_m, chr = 1){
  which(ann_m[,2] == chr)
}

getGeneStats <- function(filename, lodThreshold = 4, chrs = 1:5){
  gene_eff <- read.table(paste("Data/gene_eff~lm/", filename, sep=""))
  probe_exp <- read.table(paste("Data/gene_data/", gsub("_G", "", filename), sep=""))
  sum_num <- NULL
  idmatrix <- vector("list",5)
  res <- NULL
  nprobes <- nrow(gene_eff)
  for(chr in chrs){
    s <- 0
    ind <- NULL
    for(p in 1:nprobes){
      if(any(gene_eff[p, getProbesOnChr(ann_m, chr)] >= lodThreshold)){
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
