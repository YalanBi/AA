#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 05-02-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("C:/Arabidopsis Arrays")

peak.detect <- function(qtlprofiles, verbose = FALSE){
  cutoff = -log10(0.05 / (716 * nrow(qtlprofiles)))
  if(verbose) cat("Starting peak detection LOD >=", cutoff, "\n")
  mmatrix <- NULL
  for(x in 1:nrow(qtlprofiles)){
    peak <- FALSE
    curmax <- 0
    curmaxindex <- 1
    marker <- 1
    maximums <- NULL
    mrow <- rep(0,ncol(qtlprofiles))
    for(ab in (qtlprofiles[x,]>cutoff | qtlprofiles[x,]<(-cutoff))){
      if(ab){
        peak <- TRUE
        if(qtlprofiles[x,marker]/abs(qtlprofiles[x,marker]) > 0){
          if(qtlprofiles[x,marker] > curmax){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- marker
          }
        }else{
          if(qtlprofiles[x,marker] < (-curmax)){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- -marker
          }
        }
        if(ncol(qtlprofiles)==marker){
          if(curmax!=0) maximums <- c(maximums,curmaxindex)
        }
      }else{
        if(curmax!=0) maximums <- c(maximums,curmaxindex)
        peak <- FALSE
        curmax <- 0
      }
      marker <- marker+1
    }
    mrow[which(qtlprofiles[x,] > cutoff)] <- 1
    mrow[which(qtlprofiles[x,] < -cutoff)] <- -1
    for(a in which(maximums>0)){
      mrow[maximums[a]] <- 2
    }
    for(b in which(maximums<0)){
      mrow[(-maximums[b])] <- -2
    }
    mmatrix <- rbind(mmatrix,mrow)
  }
  rownames(mmatrix) <- rownames(qtlprofiles)
  colnames(mmatrix) <- colnames(qtlprofiles)
  #mmatrix
  write.table(mmatrix, file="Data/genes_by_chromosomes_norm_hf_cor_peak.txt")#change output filename
}



#QTL probes
qtlprofiles <- NULL
for(filename in dir("Data/")[grepl("by", dir("Data/")) & grepl("QTL", dir("Data/"))]){
  qtlprofiles <- rbind(qtlprofiles, read.table(paste("Data/", filename, sep=""), row.names=1, header=T))
}
peak.detect(qtlprofiles, T)

#good probes
qtlprofiles <- NULL
for(filename in dir("Data/")[grepl("summarized", dir("Data/")) & grepl("QTL", dir("Data/"))]){
  qtlprofiles <- rbind(qtlprofiles, read.table(paste("Data/", filename, sep=""), row.names=1, header=T))
}
peak.detect(qtlprofiles, T)
