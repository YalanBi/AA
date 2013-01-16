#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 16-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

classifyProbes <- function(filename, expmean = 5){
  res <- getGeneStats(filename, lodThreshold = 4, chrs = 1:5)
  allprobes <- 1:res$nprobes
  qtlprobes <- unique(unlist(res$ind))
  expprobes <- which(res$means >= expmean)
  goodprobes <- unique(c(qtlprobes,expprobes))
  badprobes <- which(!allprobes %in% goodprobes)
  cat(filename, badprobes,"\n")
}