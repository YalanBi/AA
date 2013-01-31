#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 29-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("C:/Arabidopsis Arrays/")

probeID <- function(np=10000){
  set.seed(1)
  ind <- sample(1:1903200, np, replace=FALSE)
  orderind <- ind[order(ind, decreasing=FALSE)]
  return(orderind)
}

getProbes <- function(np=10000){
  st <- proc.time()
  id <- probeID(np)
  cnt <- 1
  for(filename in dir("Yalan map/")[grepl("normalized",dir("Yalan map/"))]){
    cat("select from", filename, "\n")
    fp <- file(paste("Yalan map/", filename, sep=""))
    open(fp)
    mline <- ""
    header <- c("ID", strsplit(readLines(fp, n=1), "\t")[[1]])
    cat(paste(header, collapse="\t"), "\n", file="Data/new10kProbes.txt", sep="")
    mline <- readLines(fp, n=1)
    while(length(mline) != 0L & cnt <= np){
      elements <- strsplit(mline,"\t")[[1]]
      names(elements) <- header
      if(elements["ID"] %in% id){
        cat(elements["ID"], "selected\n")
        cat(mline,"\n", file="Data/new10kProbes.txt", sep="", append=TRUE)
        cnt <- cnt + 1
      }
      mline <- readLines(fp, n=1)
    }
    close(fp)
  }
  et <- proc.time()
  cat("new probes selection finishes in", (et-st)[3], "sec.\n")
}