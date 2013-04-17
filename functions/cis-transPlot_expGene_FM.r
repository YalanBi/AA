#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 15-04-2013
# first written: 15-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")

#load map file---markers' Morgan position
mpos <- read.table("refined map/map.txt")

########## for x, use Morgan ##########
totlength <- 0
lengthsx <- 0
msumlength <- NULL
for(x in 1:5){
  #msumlength <- each marker's Morgan position
  msumlength <- c(msumlength, mpos[which(mpos[,1]==x),2]+totlength)
  #totlength <- total length until current chr(Morgan)
  totlength = totlength + max(mpos[which(mpos[,1]==x),2])
  #lengthsx <- c(0, the max position of each chr separately(Morgan))
  lengthsx <- c(lengthsx, totlength)
}

########## for y, use bp ##########
#the bp length of each chr
chr1 <- 34964571
chr2 <- 22037565
chr3 <- 25499034
chr4 <- 20862711
chr5 <- 31270811
#lengthsy <- c(0, the max position of each chr separately(bp))
lengthsy <- c(0, chr1, chr1+chr2, chr1+chr2+chr3, chr1+chr2+chr3+chr4, chr1+chr2+chr3+chr4+chr5)



#qtl's
plot(c(min(lengthsx), max(lengthsx)), c(min(lengthsy), max(lengthsy)), t='n', xlab="QTL position (cM)", ylab="Gene position (bp)")
abline(h=lengthsy)
abline(v=lengthsx)
y <- 0
for(chr in 1:5){
  qtlgood <- read.table(paste("Data/genes_by_chromosomes", chr, "_norm_hf_cor_QTL.txt", sep=""), row.names=1, header=TRUE)
  location <- paste("Data/chr", chr, "_norm_hf_cor/", sep="")
  dir(location)[grepl(".png", dir(location))]
  for(fn_s in rownames(qtlgood)){
    exp <- read.table(paste("Data/chr", chr, "_norm_hf_cor/", gsub("_[12345]", "", gsub(" Response ", "", fn_s)), ".txt", sep=""), row.names=1, header=TRUE)
    y <- rep(mean(exp[ ,"bp"]), 716) + lengthsy[chr]
    colz <- as.numeric(qtlgood[fn_s,])
    col2 <- colz
    col2[colz < 5] <- rgb(1, 1, 1,0)
    col2[colz >= 5 & colz < 10] <- rgb(0.5, 0.5, 0.5 , 0.5)
    col2[colz >= 10] <- "black"
    points(x = msumlength, y = y, col = col2, pch = 20, cex = 0.5)
    cat(fn_s, "\n")
  }
}
