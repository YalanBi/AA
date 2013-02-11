


setwd("C:/Arabidopsis Arrays")

total <- 0
for(filename in dir("Data/")[grepl("by",dir("Data/")) & !grepl("QTL",dir("Data/"))]){
  chr <- which(filename == dir("Data/")[grepl("by",dir("Data/")) & !grepl("QTL",dir("Data/"))])
  aa <- read.table(paste("Data/", filename, sep=""), row.names=1)
  total <- total + nrow(aa))
}



nprobes <- nrow(aa)
plot(c(0.5, nprobes+0.5), c(round(min(aa)-0.5), round(max(aa)+0.5)))
for(p in 1:nprobes){
  points(rep(p,148), aa[p,], t='p', pch=20, cex=0.7)
}

plot(c(1, 716), c(round(min(qtl)-0.5), round(max(qtl)+0.5)))
for(p in 1:nprobes){
  points(seq(1,716), qtl[p,], t='p', pch=20, cex=0.7)
}
which.max(colMeans(qtl))
#X495326 
#    172

nn <- as.numeric(rownames(aa))

geno <- read.table("refined map/genotypes.txt",sep="\t",row.names=1,header=TRUE)
plot(c(min(nn[(nprobes-15):nprobes]), max(nn)), c(round(min(aa)-0.5), round(max(aa)+0.5)))
for(p in (nprobes-15):nprobes){
  if(-log10(anova(lm(as.numeric(aa[p,]) ~geno[,166]))[[5]][1]) > 2){
    points(rep(as.numeric(rownames(aa))[p],length(aa[p,]))+(geno[,166]/2), aa[p,], t='p', col=geno[,166], pch=20, cex=0.7)
  }
}

