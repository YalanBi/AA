#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 02-01-2013
# first written: 02-10-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("D:/Arabidopsis Arrays")

#skipping exons
seMatrix <- NULL
for(chr in 1:5){
  seMatrix <- rbind(seMatrix, read.table(paste0("Data/AS/skippingExon_chr", chr, "_wt_p3.txt"), row.names=NULL))#no NAs
}

topSEmatrix <- NULL
for(env in 1:4){
  #get top 50 SE genes: sort -log10 and take the max 50
  topSEmatrix <- cbind(topSEmatrix, paste0(seMatrix[order(seMatrix[ ,env+2], decreasing=TRUE)[1:50], 1], "_", seMatrix[order(seMatrix[ ,env+2], decreasing=TRUE)[1:50], 2]))
}
colnames(topSEmatrix) <- colnames(seMatrix)[3:6]
write.table(topSEmatrix, file="Data/AS/top50_skippingExon.txt", sep="\t")



#intron retention
riMatrix <- NULL
for(chr in 1:5){
  riMatrix <- rbind(riMatrix, read.table(paste0("Data/AS/retainedIntron_chr", chr, "_wt_p2-NA_DA.txt"), row.names=1))#NAs included
}

topRImatrix <- NULL
for(env in 1:4){
  #get top 50 RI genes: sort p-values and take the max 50, NAs are at the end of the ordered list
  topRImatrix <- cbind(topRImatrix, rownames(riMatrix)[order(riMatrix[ ,env], decreasing=TRUE)[1:50]])
}
colnames(topRImatrix) <- colnames(riMatrix)
write.table(topRImatrix, file="Data/AS/top50_retainedIntron.txt", sep="\t")



#5'&3' as site
ssMatrix <- NULL
for(chr in 1:5){
  ssMatrix <- rbind(ssMatrix, read.table(paste0("Data/AS/splicing3&5_chr", chr, "_wt_p2.txt"), row.names=1))#no NAs
}

top3Smatrix <- NULL
for(env in 1:4){
  #get top 50 3'as genes: sort p-values and take the min 50
  top3Smatrix <- cbind(top3Smatrix, rownames(ssMatrix)[order(ssMatrix[ ,env])[1:50]])
}
colnames(top3Smatrix) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
write.table(top3Smatrix, file="Data/AS/top50_splicing3.txt", sep="\t")

top5Smatrix <- NULL
for(env in 1:4){
  #get top 50 5'as genes: sort p-values and take the min 50
  top5Smatrix <- cbind(top5Smatrix, rownames(ssMatrix)[order(ssMatrix[ ,env+4])[1:50]])
}
colnames(top5Smatrix) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
write.table(top5Smatrix, file="Data/AS/top50_splicing5.txt", sep="\t")
