#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 12-04-2013
# first written: 12-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
ann_m <- read.table("refined map/map.txt")

getProbesOnChr <- function(ann_m, chr = 1){
  which(ann_m[,1] == chr)
}

#make a list of gene having at leasr 4 probes showing sig QTL(>=lodthreshould=4) on at least 1 chr
getQTLGenes <- function(chromosome = 1, empty = TRUE){  #TODO: Per chromosome we need a new file !!!
  location <- paste("D:/Arabidopsis Arrays/Data/chr",chromosome,"_norm_hf_cor/",sep="")
  
  #if file exists and we want to make it empty, empty=TRUE then empty existing file
  if(empty) cat("", file=paste("Data/genes_by_chromosomes", chromosome, "_norm_hf_cor.txt", sep=""))
  
  #load classification file, containing good probes, bad probes and so on
  load(paste(location,"Classification_chr",chromosome,"_norm_hf_cor.Rdata",sep=""))
  
  #load genomap for marker
  ann_m <- read.table("D:/Arabidopsis Arrays/refined map/map.txt")
  
  
  for(filename in dir(location)[grepl("_QTL",dir(location))]){
    qtl_data <- read.table(paste(location, filename, sep=""))
    fname <- gsub("_QTL.txt", ".txt", filename)
    exp_data <- read.table(paste(location, fname, sep=""), header=T, row.names=1)
    class <- res[[filename]]
    if(length(which(class$goodP %in% class$introP)) > 0){
      good <- class$goodP[- which(class$goodP %in% class$introP)]
    }else{
      good <- class$goodP
    }
    cnt <- 0
    for(chr in 1:5){
      marker_ids <- getProbesOnChr(ann_m, chr)
      probesa4 <- good[which(apply(qtl_data[good,marker_ids],1,max) >= 4)] #lodThreshold=4
      if(length(probesa4) >= 4){
        #Create a new gene !
        newgenename <- gsub(".txt", paste("_",chr,sep=""), fname)
        expValues <- colMeans(exp_data[probesa4, 17:ncol(exp_data)])
        cat(newgenename, expValues, "\n", file=paste("Data/genes_by_chromosomes", chromosome, "_norm_hf_cor.txt", sep=""), append=TRUE)
        cnt <- cnt+1
      }
    }
    cat("Finished with", fname, "created:", cnt, "Genes\n")
  }
}
createGenes(chromosome = 1)
