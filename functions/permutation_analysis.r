#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 15-04-2013
# first written: 15-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")


#*************************************************************** int part ***************************************************************#
#load all int permutation results into a list
int <- NULL
for(j in 1:70){
  int <- c(int, read.table(paste0("Permutation/perms_ints_", j, ".txt"))[,2])
}

#plot all results
plot(sort(unlist(int)))

#save this list into a file
cat(int, file="permutation_int.txt", sep="\n")

#make a table of FDR and corresponding int
sort(unlist(int))[length(unlist(int))*.95]
aa <- NULL
for(r in c(0.9, 0.95, 0.99)){
  aa <- rbind(aa, c(r, sort(unlist(int))[length(unlist(int))*r]))
}
colnames(aa) <- c("FDR"," Int")
write.table(aa, file="permutation_int_FDR.txt", row.names=F)

#plot histogram of int permutation
png(file="histogram of Int permutation.png", width=480, height=480, bg="white")
hist(int, breaks=20, border=grey(0.5), col="grey", main="Histogram of Permutation for Interaction QTL", cex.main=1.5, xlab="Interaction QTL", cex.lab=1)
abline(v=11.6125174131782, col="red")
axis(1, at=11.61, labels=11.6, col="red", col.axis="red", cex.axis=1.5)
dev.off()


#*************************************************************** qtl part ***************************************************************#
#load all qtl permutation results into a list
qtl <- NULL
for(j in 1:70){
  qtl <- c(qtl, read.table(paste0("Permutation/perms_qtls_", j, ".txt"))[,2])
}

#plot all results
plot(sort(unlist(qtl)))

#save this list into a file
cat(qtl, file="permutation_qtl.txt", sep="\n")

#make a table of FDR and corresponding qtl
sort(unlist(qtl))[length(unlist(qtl))*.95]
bb <- NULL
for(r in c(0.9, 0.95, 0.99)){
  bb <- rbind(bb, c(r, sort(unlist(qtl))[length(unlist(qtl))*r]))
}
colnames(bb) <- c("FDR"," QTL")
write.table(bb, file="permutation_qtl_FDR.txt", row.names=F)

#plot histogram of qtl permutation
png(file="histogram of QTL permutation.png", width=480, height=480, bg="white")
hist(qtl, breaks=20, border=grey(0.5), col="grey", main="Histogram of Permutation for eQTL", cex.main=1.5, xlab="eQTL", cex.lab=1)
abline(v=7.95960219666518, col="red")
axis(1, at=7.96, labels=8.0, col="red", col.axis="red", cex.axis=1.5)
dev.off()



#************************************************* count nQTL/nInt higher than threshold *************************************************#
#first, QTL counting
#qtlFDR <- unlist(read.table("Permutation/permutation_qtl_FDR.txt", row.names=1, header=T))
qtlFDR <- c(7.6, 8.0, 8.8)
cat("QTL counting for 3 FDRs starts!\n")

nSample <- c(0, 0, 0)
nExon <- c(0, 0, 0)
nGene <- c(0, 0, 0)
for(chr in 1:5){
  cat("chr", chr, "starts counting...\n")
  
  countm <- read.table(paste0("Data/fullModeMapping/expGenes_chr", chr, "_FMD_QTL.txt"), row.names=1, header=T)
  
  for(fdr in 1:3){
    nSample[fdr] <- nSample[fdr] + length(which(abs(countm) >= qtlFDR[fdr], arr.ind=F))
    cat("higher than", qtlFDR[fdr], ":", length(which(abs(countm) >= qtlFDR[fdr], arr.ind=F)), "RILs\n")
    nExon[fdr] <- nExon[fdr] + length(unique(rownames(which(abs(countm) >= qtlFDR[fdr], arr.ind=T))))
    cat("higher than", qtlFDR[fdr], ":", length(unique(rownames(which(abs(countm) >= qtlFDR[fdr], arr.ind=T)))), "exons\n")
    nGene[fdr] <- nGene[fdr] + length(unique(apply(as.matrix(unique(rownames(which(abs(countm) >= qtlFDR[fdr], arr.ind=T)))), 1, function(x) strsplit(x , "_")[[1]][1])))
    cat("higher than", qtlFDR[fdr], ":", length(unique(apply(as.matrix(unique(rownames(which(abs(countm) >= qtlFDR[fdr], arr.ind=T)))), 1, function(x) strsplit(x , "_")[[1]][1]))), "genes\n")
  }
  cat("chr", chr, "counting finished\n")
}
res <- cbind(qtlFDR, nSample, nExon, nGene)
colnames(res) <- c("-log10(QTL)", "nSample", "nExon", "nGene")


#next, Int counting
#intFDR <- unlist(read.table("Permutation/permutation_int_FDR.txt", row.names=1, header=T))
intFDR <- c(10.4, 11.6, 13.3)
cat("Int counting for 3 FDRs starts!\n")

nSample <- c(0, 0, 0)
nExon <- c(0, 0, 0)
nGene <- c(0, 0, 0)
for(chr in 1:5){
  cat("chr", chr, "starts counting...\n")
  
  countm <- read.table(paste0("Data/fullModeMapping/expGenes_chr", chr, "_FMD_Int.txt"), row.names=1, header=T)
  
  for(fdr in 1:3){
    nSample[fdr] <- nSample[fdr] + length(which(abs(countm) >= intFDR[fdr], arr.ind=F))
    cat("higher than", intFDR[fdr], ":", length(which(abs(countm) >= intFDR[fdr], arr.ind=F)), "RILs\n")
    nExon[fdr] <- nExon[fdr] + length(unique(rownames(which(abs(countm) >= intFDR[fdr], arr.ind=T))))
    cat("higher than", intFDR[fdr], ":", length(unique(rownames(which(abs(countm) >= intFDR[fdr], arr.ind=T)))), "exons\n")
    nGene[fdr] <- nGene[fdr] + length(unique(apply(as.matrix(unique(rownames(which(abs(countm) >= intFDR[fdr], arr.ind=T)))), 1, function(x) strsplit(x , "_")[[1]][1])))
    cat("higher than", intFDR[fdr], ":", length(unique(apply(as.matrix(unique(rownames(which(abs(countm) >= intFDR[fdr], arr.ind=T)))), 1, function(x) strsplit(x , "_")[[1]][1]))), "genes\n")
  }
  cat("chr", chr, "counting finished\n")
}
res <- cbind(intFDR, nSample, nExon, nGene)
colnames(res) <- c("-log10(Int)", "nSample", "nExon", "nGene")
