#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 24-01-2013
# first written: 24-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

checkTUmeans <- function(filename, location = "C:/Arabidopsis Arrays/Data/chr1/"){
  rawexp <- read.table(paste(location, filename, sep=""), header=TRUE, row.names=1)
  newexp <- rawexp[,16:163]
  s <- 1
  tu <- NULL
  tumeans <- NULL
  for(p in 1:nrow(rawexp)){
    if(rawexp[p,"tu"] != rawexp[s,"tu"]){
      tu <- c(as.character(rawexp[s,"tu"]),s,p-1,mean(as.matrix(newexp[s:(p-1),])))
      tumeans <- rbind(tumeans, tu)
      s <- p
    }
  }
  tu <- c(as.character(rawexp[s,"tu"]),s,p,mean(as.matrix(newexp[s:p,])))
  tumeans <- rbind(tumeans, tu)
  return(tumeans)
}

checkTUmeans("AT1G01010.txt")
