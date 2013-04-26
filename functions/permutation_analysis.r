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
png(file="histogram of Int permutation.png", width=1024, height=1024, bg="white")
hist(int, breaks=100, col="darkseagreen2", border="darkolivegreen", main="Histogram of Int Permutation", cex.main=2, xlab="Int", cex.lab=1.5)
abline(v=11.6125174131782, col="darkolivegreen")
text(x=12.15, y=35, labels="FDR=95%", col="darkolivegreen", cex=1.2)
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
png(file="histogram of QTL permutation.png", width=1024, height=1024, bg="white")
hist(qtl, breaks=100, col="lightpink2", border="maroon", main="Histogram of QTL Permutation", cex.main=2, xlab="QTL", cex.lab=1.5)
abline(v=7.95960219666518, col="maroon")
text(x=9, y=125, labels="FDR=95%", col="maroon", cex=1.2)
dev.off()
