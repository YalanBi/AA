#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 29-01-2013
# first written: 24-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

load("C:/Arabidopsis Arrays/Data/OLD_unnormalized/new genes/mapNewQTLgenes.Rdata")
plotCheckTransP <- function(chrQTL){
  plot(c(1, 720), c(1, 70))
  for(chr in 1:5){
    for(x in 1:nrow(chrQTL[[paste("chr", chr, sep="")]])){
      points(1:716, chrQTL[[paste("chr", chr, sep="")]][x,], pch=20, col=chr)
    }
    axis(1, at=which.max(chrQTL$chrMeansQTL[chr,]), labels=which.max(chrQTL$chrMeansQTL[chr,]))
  }
}
plotCheckTransP(chrQTL)

plotCheckTransL <- function(chrQTL$chrMeansQTL){
  plot(c(1,720), c(1,max(chrQTL$chrMeansQTL)))
  for(chr in 1:5){
    points(chrQTL$chrMeansQTL[chr,], t='l', col=chr)
    axis(1, at=which.max(chrQTL$chrMeansQTL[chr,]), labels=which.max(chrQTL$chrMeansQTL[chr,]))
  }
}
plotCheckTransL(chrQTL)
