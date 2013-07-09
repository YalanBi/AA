#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 02-07-2013
# first written: 21-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#compare wilcoxon vs ANOVA using
setwd("D:/Arabidopsis Arrays")
old <- NULL
new <- NULL
for(chr in 1:5){
  old <- rbind(old, read.table(paste0("Data/AS/SE_chr", chr, "_wt_less.txt"), row.names=NULL)[-2])
  new <- rbind(new, read.table(paste0("Data/AS/SE_chr", chr, "_ANOVA_less.txt"), row.names=NULL)[-2])
}
png(filename="Wilcox_vs_ANOVA.png", width=1024, height=1024, bg="white")
plot(c(-1, max(rbind(old[,2:5], new[,2:5]))), c(-1, max(rbind(old[,2:5], new[,2:5]))), type='n', xlab="Wilcoxon", ylab="ANOVA")
for(e in 1:nrow(old)){
  for(i in 2:5){
    points(x=old[e,i], y=new[e,i], pch=20, cex=0.5)
  }
}
dev.off()
