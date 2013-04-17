#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 15-04-2013
# first written: 15-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
chrEdges <- c(1, 162, 263, 403, 527)


#*********************************************************** peak detection part ***********************************************************#
peakDetect <- function(data, chrEdges, cutoff = 4){
  above   <- as.numeric(data > cutoff)
  below   <- -as.numeric((data < (-cutoff)))
  pattern <- above+below

  inPeek  <- 0  # Are we in a peak ?
  top     <- 0  # What was our previous TOP score ?
  for(x in 1:length(pattern)){
    if(x %in% chrEdges){
      if(inPeek == -1) pattern[topx] <- -2
      if(inPeek == 1) pattern[topx] <- 2
      inPeek <- 0
    }
    cat(x, pattern[x], data[x], inPeek,"\n")
    if(pattern[x] == 1){
      if(inPeek == 0){
        inPeek <- 1
        top  <- data[x]
        topx <- x
      }else if(inPeek == 1){ # We're in a peak already check this one is the 'top'
        if(top < data[x]){ top <- data[x] ; topx <- x }
        if(top > data[x]){ pattern[topx] <- 2; top <- data[x]; } # Going down, might be a double peek
      }else{ stop(paste("Positive peak to negative peak at", x)) }
    }else if(pattern[x] == -1){
      if(inPeek == 0){
        inPeek <- (-1)
        top  <- data[x]
        topx <- x
      }else if(inPeek == -1){ # We're in a peak already check this one is the 'top'
        if(top > data[x]){ top <- data[x]; topx <- x } # Its a peak
        if(top < data[x]){ pattern[topx] <- -2; top <- data[x]; } # Going down, might be a double peek
      }else{ stop(paste("Negative peak to positive peak at", x)) }
    }else if(pattern[x] == 0){
      if(inPeek == 1){ inPeek = 0; pattern[topx] <- 2 }
      if(inPeek == -1){ inPeek = 0; pattern[topx] <- -2 }
    }
  }
  pattern
}


#*************************************************************** plot part ***************************************************************#
pdmatrix <- NULL
for(p in 1:nrow(qtl)){
  pdmatrix <- rbind(pdmatrix, peakDetect(data = unlist(qtl[p, ]), chrEdges, cutoff = 4))
}
pdmatrix

plot(c(1,39), c(0,610))
for(p in 1:39){
  points(x=rep(p, 716), y=msumlength, col=pdmatrix[p, ], pch=20, cex=0.5)
}


#load map file---markers' Morgan position
mpos <- read.table("refined map/map.txt")

########## for x, use Morgan ##########
totlength <- 0
lengthsx <- 0
msumlength <- NULL
for(x in 1:5){
  #msumlength <- each marker's Morgan position
  msumlength <- c(msumlength, mpos[which(mpos[,1]==x),2] + totlength + 25*(x-1))
  #totlength <- total length until current chr(Morgan)
  totlength = totlength + max(mpos[which(mpos[,1]==x),2])
  #lengthsx <- c(0, the max position of each chr separately(Morgan))
  lengthsx <- c(lengthsx, totlength)
}



for(p in 1:39){
  for(m in 1:716){
    if(pdmatrix[p, m] == 1){
      points(x=p, y=msumlength[m], col="orange", pch=20, cex=0.5)
    }
  }
  for(m in 1:716){
    if(pdmatrix[p, m] == 2){
      points(x=p, y=msumlength[m], col="darkorange3", pch=20, cex=0.5)
    }
  }
}