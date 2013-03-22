#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 22-03-2013
# first written: 22-03-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("D:/Arabidopsis Arrays")

createGenesByTU <- function(chromosome = 1, empty = TRUE){#empty the existed contents
  location <- paste("D:/Arabidopsis Arrays/Data/chr", chromosome, "_norm_hf_cor/", sep="")
  output <- paste("Data/genesByTU_chr", chromosome, "_norm_hf_cor.txt", sep="")
  if(empty) cat("", file=output)
  load(paste(location,"Classification_chr",chromosome,"_norm_hf_cor.Rdata",sep=""))
  
  header_idx <- 1
  for(filename in dir(location)[grepl("_QTL",dir(location))]){
    fname <- gsub("_QTL.txt",".txt", filename)
    exp_data <- read.table(paste(location, fname, sep=""),header=T,row.names=1)
    class <- res[[filename]]#res[["AT1G01010_QTL.txt"]]
    
    if(header_idx == 1){
      cat(c("ID", colnames(exp_data)[17:ncol(exp_data)], "\n"), file =output)
      header_idx <- 2
    }
    
    if(length(which(class$goodP %in% class$introP)) > 0){
      good <- class$goodP[- which(class$goodP %in% class$introP)]
    }else{
      good <- class$goodP
    }
    
    if(length(good) > 0){
      for(tu in gsub("tu", "", unique(exp_data[good,"tu"]))){
        expValues <- colMeans(exp_data[good[which(exp_data[good,"tu"]==paste("tu", tu, sep=""))],17:ncol(exp_data)])
        newgenename <- paste(gsub(".txt", "", fname), "_tu", tu, sep="")
        cat(newgenename, expValues, "\n", file=output, append=TRUE)
      }
      cat("Finished with", fname, "\n")
    }else{
      cat(fname, "no good probes\n")
    }
  }
}
createGenesByTU(chromosome = 1)
