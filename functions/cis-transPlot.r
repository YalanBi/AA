#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 29-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("C:/Arabidopsis Arrays")
mpos <- read.table("refined map/map.txt")
chr1Length <- max(mpos[which(mpos[,1]==1),])
chr2Length <- max(mpos[which(mpos[,1]==2),])
chr3Length <- max(mpos[which(mpos[,1]==3),])
chr4Length <- max(mpos[which(mpos[,1]==4),])
chr5Length <- max(mpos[which(mpos[,1]==5),])
x <- c(mpos[which(mpos[,1]==1),2], mpos[which(mpos[,1]==2),2]+chr1Length, mpos[which(mpos[,1]==3),2]+chr1Length+chr2Length, mpos[which(mpos[,1]==4),2]+chr1Length+chr2Length+chr3Length, mpos[which(mpos[,1]==5),2]+chr1Length+chr2Length+chr3Length+chr4Length)

aa <- read.table("Data/chr1_normalized/AT1G01010.txt", row.names=1, header=TRUE)
aap <- aa[,17:ncol(aa)]
chraa <- read.table("Data/genes_summarized_1_normalized.txt", row.names=1, header=FALSE)
for(chr in 1:5){
  qtlchr <- read.table(paste("Data/genes_summarized_", chr, "_normalized_QTL.txt", sep=""), row.names=1, header=TRUE)
  plot(c(0, max(x)), )
  for(fn in rownames(qtlchr)){
    exp <- read.table(paste("Data/chr1_normalized/", fn, ".txt", sep=""))
    pForCol <- round(qtlchr[fn, ]-2.5)
    pForCol[pForCol > 8] <- 8
    pForCol[pForCol < 1] <- 1
    points(x, rep(means(exp[ ,"bp"]), 716), col=OrRd[unlist(pForCol)], pch=15)
    
    }
    line()
  }
}

OrRd <- brewer.pal(8,"YlOrRd")





expchr <- read.table(paste("Data/genes_summarized_", chr, "_normalized.txt", sep=""), row.names=1, header=FALSE)

