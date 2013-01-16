#P5-distribution of rarios
setwd("X:/Data/Desktop/Arabidopsis Arrays")
ratio <- read.table("Data/ratio.txt")
#make a plot
png(filename="Data/ratio_distribution.png", width = 600, height = 600)
plot(c(1,5), c(0,1), xlab="chr", ylab="ratio", cex.axis=0.75, cex.lab=1, las=1, mgp=c(1.75,0.75,0), t="n")
for(chr in 1:ncol(ratio)){
  points(rep(chr, (nrow(ratio)-1))+((runif(nrow(ratio)-1)-0.5)/15), ratio[2:nrow(ratio),chr], pch=20, cex=0.5, col=rgb(0.4,0,0.6,0.25))
}
dev.off()


#make a boxplot + points
boxplot(ratio,boxwex = 0.15)
for(chr in 1:ncol(ratio)){
  points(rep(chr, (nrow(ratio)-1))+((runif(nrow(ratio)-1)-0.5)/15), ratio[2:nrow(ratio),chr], pch=20, cex=0.5, col=rgb(0.4,0,0.6,0.2))
}

#make a .txt file
perc <- NULL
for(chr in 1:ncol(ratio)){
  a <- 0
  b <- 0
  c <- 0
  d <- 0
  e <- 0
  f <- 0
  for(g in 1:nrow(ratio)){
    if(ratio[g,chr]>0 & ratio[g,chr]<=0.1){
      a <- a+1
    } else if(ratio[g,chr]>0.1 & ratio[g,chr]<=0.2){
      b <- b+1
    } else if(ratio[g,chr]>0.2 & ratio[g,chr]<=0.3){
      c <- c+1
    } else if(ratio[g,chr]>0.3 & ratio[g,chr]<=0.4){
      d <- d+1
    } else if(ratio[g,chr]>0.4 & ratio[g,chr]<=0.5){
      e <- e+1
    } else{
      f <- f+1
    }
  }
  perc <- cbind(perc,c(a,b,c,d,e,f)/nrow(ratio))
}
rownames(perc) <- c("(0,0.1]","(0.1,0.2]","(0.2,0.3]","(0.3,0.4]","(0.4,0.5]","(0.5,1]")
colnames(perc) <- paste("chr", 1:5, sep="")
write.table(perc, file="Data/ratio_distribution.txt")