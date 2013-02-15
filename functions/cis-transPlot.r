#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 4-02-2013
# first written: 1-02-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("C:/Arabidopsis Arrays")
mpos <- read.table("refined map/map.txt")

totlength <- 0
lengthsx <- 0
for(x in 1:5){
  totlength = totlength + max(mpos[which(mpos[,1]==x),2])
  lengthsx <- c(lengthsx,totlength)
}


chr1 <- 34964571
chr2 <- 22037565
chr3 <- 25499034
chr4 <- 20862711
chr5 <- 31270811
lengthsy <- c(0, chr1, chr1+chr2, chr1+chr2+chr3, chr1+chr2+chr3+chr4, chr1+chr2+chr3+chr4+chr5)


chr1Length <- max(mpos[which(mpos[,1]==1),])
chr2Length <- max(mpos[which(mpos[,1]==2),])
chr3Length <- max(mpos[which(mpos[,1]==3),])
chr4Length <- max(mpos[which(mpos[,1]==4),])
chr5Length <- max(mpos[which(mpos[,1]==5),])
msumlength <- c(mpos[which(mpos[,1]==1),2], mpos[which(mpos[,1]==2),2]+chr1Length, mpos[which(mpos[,1]==3),2]+chr1Length+chr2Length, mpos[which(mpos[,1]==4),2]+chr1Length+chr2Length+chr3Length, mpos[which(mpos[,1]==5),2]+chr1Length+chr2Length+chr3Length+chr4Length)

#good's
plot(c(min(lengthsx),max(lengthsx)),c(min(lengthsy),max(lengthsy)),t='n',xlab="QTL position (cM)", ylab="Gene position (bp)")
abline(h=lengthsy)
abline(v=lengthsx)
y <- 0
for(chr in 1:5){
  qtlchr <- read.table(paste("Data/genes_summarized_", chr, "_norm_hf_cor_QTL.txt", sep=""), row.names=1, header=TRUE)
  for(fn in rownames(qtlchr)){
    exp <- read.table(paste("Data/chr", chr, "_norm_hf_cor/", gsub(" Response ", "", fn), ".txt", sep=""), row.names=1, header=TRUE)
    y <- rep(mean(exp[ ,"bp"]), 716) + lengthsy[chr]
    colz <- as.numeric(qtlchr[fn,])
    col2 <- colz
    col2[colz < 5] <- rgb(1, 1, 1,0)
    col2[colz >= 5 & colz < 10] <- rgb(0.5, 0.5, 0.5 , 0.5)
    col2[colz >= 10] <- "black"
    points(x=msumlength, y=y, col=col2, pch=20,cex=0.5)
    cat(fn, "\n")
  }
}


#qtl's
plot(c(min(lengthsx),max(lengthsx)),c(min(lengthsy),max(lengthsy)),t='n',xlab="QTL position (cM)", ylab="Gene position (bp)")
abline(h=lengthsy)
abline(v=lengthsx)
y <- 0
for(chr in 1:5){
  qtlgood <- read.table(paste("Data/genes_by_chromosomes", chr, "_norm_hf_cor_QTL.txt", sep=""), row.names=1, header=TRUE)
  location <- paste("Data/chr", chr, "_norm_hf_cor/", sep="")
  dir(location)[grepl(".png", dir(location))]
  for(fn_s in rownames(qtlgood)){
    exp <- read.table(paste("Data/chr", chr, "_norm_hf_cor/", gsub("_[12345]","", gsub(" Response ", "", fn_s)), ".txt", sep=""), row.names=1, header=TRUE)
    y <- rep(mean(exp[ ,"bp"]), 716) + lengthsy[chr]
    colz <- as.numeric(qtlgood[fn_s,])
    col2 <- colz
    col2[colz < 5] <- rgb(1, 1, 1,0)
    col2[colz >= 5 & colz < 10] <- rgb(0.5, 0.5, 0.5 , 0.5)
    col2[colz >= 10] <- "black"
    points(x=msumlength, y=y, col=col2, pch=20,cex=0.5)
    cat(fn_s, "\n")
  }
}











chrslength_x <- NULL
lastchrlength <- 0
x <- NULL
chrslength_y <- NULL
ymax <- 0
for(chr in 1:5){
  chrslength_x <- c(chrslength_x, max(mpos[which(mpos[,1]==chr),]))
  mp <- mpos[which(mpos[,1]==chr),2] + lastchrlength
  x <- c(x, mp)
  lastchrlength <- lastchrlength + chrslength_x[chr]
  location <- paste("Data/chr", chr, "_normalized/", sep="")
  aa <- read.table(paste(location, max(dir(location)[grepl(".txt",dir(location)) & !grepl("_QTL",dir(location))]), sep=""), row.names=1, header=TRUE)
  chrslength_y <- c(chrslength_y, mean(aa[ ,"bp"]))
  ymax <- ymax + chrslength_y
}
plot(c(0, max(x)), c(0, y), col='white')
for(chr in 1:5){
  st <- proc.time()
  qtlchr <- read.table(paste("Data/genes_summarized_", chr, "_normalized_QTL.txt", sep=""), row.names=1, header=TRUE)
  
  for(fn in rownames(qtlchr)){
    exp <- read.table(paste("Data/chr", chr, "_normalized/", gsub(" Response ", "", fn), ".txt", sep=""), row.names=1, header=TRUE)
    pForCol <- round(qtlchr[fn, ]-1.5)
    pForCol[pForCol > 9] <- 9
    pForCol[pForCol < 2] <- 2
    points(x, rep(mean(exp[ ,"bp"]), 716), col=OrRd[unlist(pForCol)], pch=15)
    cat(fn, "\n")
  }
  line()
  et <- proc.time()
  cat("chr", chr, "finished in", (et-st)[3], "secs\n")
}
