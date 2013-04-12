#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 09-04-2013
# first written: 22-03-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("D:/Arabidopsis Arrays")

#calculate the persentage of introns(number of genes containing retented introns divided by total number of genes of one chromosome, each chr separately)
IntronRetentionPost <- function(threshould = 7){
  ratiomatrix <- NULL
  
  for(chr in 1:5){
    #load intron retention analysis file
    aa <- read.table(paste0("Data/intronRetention_chr", chr, ".txt"), row.names=1, header=T)
    #calculate gene number of this chr
    UniqueGenes <- unique(as.character(unlist(lapply(strsplit(rownames(aa),"_"),"[[",1))))
    
    ratio <- NULL
    
    #each chr separately
    for(env in 1:4){
      #retained <- the names of all retented introns of this chr
      retained <- rownames(aa[which(aa[ ,(env*4-1)] < threshould),])
      #uniqueRetained <- the names of genes containing retented introns of this chr
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






#calculate the persentage of introns(number of genes containing retented introns divided by total number of genes of all 5 chromosome)
IntronRetentionPostAll <- function(threshould = 7){
  #reline <- number of genes containing retented introns for all 5 chr, under 4 Env separately
  retainedline <- c(0, 0, 0, 0)
  #ngenes <- the sum of genes of 5 chr
  ngenes <- 0
  
  for(chr in 1:5){
    #load intron retention analysis file
    aa <- read.table(paste0("Data/intronRetention_chr", chr, ".txt"), row.names=1, header=T)
    #uniqueRetained <- the names of genes containing retented introns of this chr
    UniqueGenes <- unique(as.character(unlist(lapply(strsplit(rownames(aa),"_"),"[[",1))))
    ngenes <- ngenes + length(UniqueGenes)
    
    #reline <- number of genes containing retented introns, for this chr, under 4 Env separately
    reline <- NULL
    for(env in 1:4){
      #for this chr, under this Env, the names of retented introns
      retained <- rownames(aa[which(aa[ ,(env*4-1)] < threshould),])
      #for this chr, under this Env, the names of genes containing retented introns
      UniqueRetained <- unique(as.character(unlist(lapply(strsplit(retained,"_"),"[[",1))))
      
      reline <- c(reline, length(UniqueRetained))
    }
    retainedline <- retainedline + reline
  }
  ratio <- retainedline/ngenes
  
  return(ratio)
}

IntronRetentionPostAll(threshould = 7)
