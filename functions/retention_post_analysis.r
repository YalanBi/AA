#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 26-03-2013
# first written: 22-03-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("D:/Arabidopsis Arrays")

IntronRetentionPost <- function(threshould = 7){
  ratiomatrix <- NULL
  
  for(chr in 1:5){
    aa <- read.table(paste0("Data/intronRetention_chr", chr, ".txt"), row.names=1, header=T)
    UniqueGenes <- unique(as.character(unlist(lapply(strsplit(rownames(aa),"_"),"[[",1))))
    ratio <- NULL
    
    for(env in 1:4){
      retained <- rownames(aa[which(aa[ ,(env*4-1)] < threshould),])
      UniqueRetained <- unique(as.character(unlist(lapply(strsplit(retained,"_"),"[[",1))))
      
      ratio <- c(ratio, length(UniqueRetained) / length(UniqueGenes))
    }
    
    ratiomatrix <- rbind(ratiomatrix, ratio)
  }
  
  rownames(ratiomatrix) <- c("chr1", "chr2", "chr3", "chr4", "chr5")
  colnames(ratiomatrix) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
  
  return(ratiomatrix)
}

retention <- IntronRetentionPost(threshould = 7)

tres <- NULL
for(a in 1:4){
  tres <- c(tres, -log10(t.test(retention[ ,a])$p.value))
}
retention <- rbind(retention, tres)


IntronRetentionPostAll <- function(threshould = 7){
  retainedline <- c(0, 0, 0, 0)
  ngenes <- 0
  for(chr in 1:5){
    aa <- read.table(paste0("Data/intronRetention_chr", chr, ".txt"), row.names=1, header=T)
    UniqueGenes <- unique(as.character(unlist(lapply(strsplit(rownames(aa),"_"),"[[",1))))
    ngenes <- ngenes + length(UniqueGenes)
    reline <- NULL
    for(env in 1:4){
      retained <- rownames(aa[which(aa[ ,(env*4-1)] < threshould),])
      UniqueRetained <- unique(as.character(unlist(lapply(strsplit(retained,"_"),"[[",1))))
      
      reline <- c(reline, length(UniqueRetained))
    }
    retainedline <- retainedline + reline
  }
  ratio <- retainedline/ngenes
  
  return(ratio)
}

IntronRetentionPostAll(threshould = 7)
