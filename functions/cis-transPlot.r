#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 10-04-2013
# first written: 01-02-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

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


#for chr1/chr5, use function "loadQTLfile" to load the QTL_t file
loadQTLfile <- function(chr = 1){
  x <- 0
  matrix <- NULL
  mnames <- NULL
  fp <- file(paste0("Data/newTU/genesByTU_chr", chr, "_norm_hf_cor_QTL_t.txt"))
  open(fp)
  line <- readLines(fp, 1)
  tnames <- strsplit(line,"\t")[[1]][-1]
  line <- readLines(fp, 1)
  x <- x+1
  
  while(length(line) > 0 && line != ""){
    data <- strsplit(line,"\t")[[1]]
    mnames <- c(mnames, data[1])
    #cat("Done", x, "\n")
    
    data <- round(as.numeric(data[-1]), d=2)
    matrix <- cbind(matrix, data)
    line <- readLines(fp, 1)
    
    x <- x+1
  }
  colnames(matrix) <- mnames
  rownames(matrix) <- tnames
  
  close(fp)
  invisible(matrix)
}


drawPoints <- function(qtlchr){
  #all TUs' names of this chr
  tunames <- rownames(qtlchr)
  #all genes' names of this chr
  genenames <- unique(as.character(unlist(lapply(strsplit(rownames(qtlchr),"_"),"[[",1))))
  for(fn in genenames){
    #load file "AT1G01010.txt"
    cat("loading", fn, "file\n")
    exp <- read.table(paste("Data/chr", chr, "_norm_hf_cor/", fn, ".txt", sep=""), row.names=1, header=TRUE)
    
    for(tu in tunames[which(grepl(fn, tunames))]){
      y <- rep(mean(exp[which(exp[,"tu"]==strsplit(tu, "_")[[1]][2]),"bp"]), 716) + lengthsy[chr]
      cat(tu, "bp is:", unique(y), "\n")
      
      colz <- as.numeric(qtlchr[tu,])
      col2 <- colz
      col2[colz < 5] <- rgb(1, 1, 1,0)
      col2[colz >= 5 & colz < 10] <- rgb(0.5, 0.5, 0.5 , 0.5)
      col2[colz >= 10] <- "black"
      points(x=msumlength, y=y, col=col2, pch=20,cex=0.5)
    }
  }
}



#geneByTU
plot(c(min(lengthsx), max(lengthsx)), c(min(lengthsy), max(lengthsy)), t='n', xlab="QTL position (cM)", ylab="Gene position (bp)")
abline(h=lengthsy)
abline(v=lengthsx)
y <- 0
for(chr in 1:5){
  qtlchr <- NULL
  if(chr == 1 | chr == 5){
    cat("loading chr", chr, "QTL_t file...\n")
    qtlchr <- loadQTLfile(chr = chr)
  }else{
    cat("loading chr", chr, "QTL file...\n")
    qtlchr <- read.table(paste("Data/newTU/genesByTU_chr", chr, "_norm_hf_cor_QTL.txt", sep=""), row.names=1, header=TRUE)
  }
  tunames <- rownames(qtlchr)
  #all genes' names of this chr
  genenames <- unique(as.character(unlist(lapply(strsplit(rownames(qtlchr),"_"),"[[",1))))
  for(fn in genenames){
    #load file "AT1G01010.txt"
    cat("loading", fn, "file\n")
    expressions <- read.table(paste("Data/chr", chr, "_norm_hf_cor/", fn, ".txt", sep=""), row.names=1, header=TRUE)
    
    for(tu in tunames[which(grepl(fn, tunames))]){
      y <- rep(mean(expressions[which(expressions[,"tu"]==strsplit(tu, "_")[[1]][2]),"bp"]), 716) + lengthsy[chr]
      cat(tu, "bp is:", unique(y), "\n")
      
      colz <- as.numeric(qtlchr[tu,])
      col2 <- colz
      col2[colz < 5] <- rgb(1, 1, 1,0)
      col2[colz >= 5 & colz < 10] <- rgb(0.5, 0.5, 0.5 , 0.5)
      col2[colz >= 10] <- "black"
      points(x=msumlength, y=y, col=col2, pch=20,cex=0.5)
    }
    rm(expressions)
    gc();gc();gc();gc();gc();gc();gc();
  }
  rm(qtlchr)
  gc();gc();gc();gc();gc();gc();gc();
}



##########################################useless######################################
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
    exp <- read.table(paste("Data/chr", chr, "_norm_hf_cor/", gsub("_[12345]", "", gsub(" Response ", "", fn_s)), ".txt", sep=""), row.names=1, header=TRUE)
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





###############################################OLD|OLD|OLD###########################################
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
