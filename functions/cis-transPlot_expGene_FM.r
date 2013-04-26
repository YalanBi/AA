#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 15-04-2013
# first written: 15-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")

#load map file---markers' Morgan position
mpos <- read.table("refined map/map_bp_remove_nogap.txt", row.names=1, header=T)


#for x, use marker bp; for y, use probe bp
totlength <- c(0, 34964571, 22037565, 25499034, 20862711, 31270811)
chrEnds <- c(161, 262, 402, 526, 716)
lpos <- NULL
for(chr in 1:5){
  #add gap between chrs, then get markers' bp position!
  mpos[which(mpos[,1] == chr),2] <- (mpos[which(mpos[,1] == chr),2] + 6000000*(chr-1))/1000
  lpos <- c(lpos, (sum(totlength[1:chr]) + 6000000*(chr-1))/1000, (sum(totlength[1:(chr + 1)]) + 6000000*(chr-1))/1000)
}
#lpos(0.00, 34964.57, 40964.57, 63002.14, 69002.14, 94501.17, 100501.17, 121363.88, 127363.88, 158634.69)


#************************************************************** plot cis-trans **************************************************************#
lodthreshold = 8
lodthreshold_black = 12

png(filename = paste0("cis-trans_exon_", lodthreshold, "+", lodthreshold_black, ".png"), width = 1024, height = 1024, bg = "white")
plot(x = c(5000, max(mpos, na.rm=T)-1000), y = c(5000, max(lpos)-5000), t='n', xlab="QTL position (kb)", ylab="Gene position (kb)", xaxt='n', yaxt='n')
axis(1, at=seq(0, 150000, 10000), labels=seq(0, 150000, 10000))
axis(2, at=seq(0, 150000, 10000), labels=seq(0, 150000, 10000))

for(r in 1:4){
  rect(lpos[r*2], -9000000, lpos[r*2+1], max(lpos)+9000000, col=rgb(159, 182, 205, 105, max = 255), border="transparent")
  rect(-9000000, lpos[r*2], max(lpos)+9000000, lpos[r*2+1], col=rgb(159, 182, 205, 105, max = 255), border="transparent")
}
#abline(h=lpos[seq(1,10,2)], v=lpos[seq(1,10,2)],col='green',lty=2)
#abline(h=lpos[seq(2,10,2)], v=lpos[seq(2,10,2)],col='red',lty=2)

for(chr in 1:5){
  st <- proc.time()
  
  cat("chr", chr, "position loading...\n")
  #expbp <- read.table(paste("Data/fullModeMapping/expGenes_chr", chr, "_1tu1probe.txt", sep=""), row.names=1, header=TRUE)[,"bp"]
  #another way to load exon location (bp)
  fbp <- file(paste("Data/fullModeMapping/expGenes_chr", chr, "_1tu1probe.txt", sep=""))
  open(fbp)
  invisible(readLines(fbp, 1))
  pline <- readLines(fbp, 1)
  expbp <- NULL
  while(length(pline) > 0 && pline != ""){
    expbp <- c(expbp, as.numeric(strsplit(pline,"\t")[[1]][2])/1000)
    pline <- readLines(fbp, 1)
  }
  close(fbp)
  et <- proc.time()
  cat("chr", chr, "got bp after", (et-st)[3], "secs!\n")
  
  
  fp <- file(paste("Data/fullModeMapping/expGenes_chr", chr, "_FM_QTL.txt", sep=""))
  open(fp)
  mline <- readLines(fp, 1)
  
  for(p in 1:length(expbp)){
    mline <- readLines(fp, 1)
    qtl <- as.numeric(strsplit(mline,"\t")[[1]][2:717])
    col2 <- qtl
    col2[qtl < lodthreshold] <- rgb(1, 1, 1,0)
    col2[qtl >= lodthreshold & qtl < lodthreshold_black] <- "gray80"
    col2[qtl >= lodthreshold_black] <- "black"
    points(x = mpos[,2], y = rep(expbp[p], 716) + (sum(totlength[1:chr]) + 6000000*(chr-1))/1000, col = col2, pch = 20, cex = 0.5)
  }
  
  close(fp)
  
  et <- proc.time()
  cat("chr", chr, "finished plotting after", (et-st)[3], "secs!\n")
}
dev.off()





#********************************************************** remove trans markers **********************************************************#
#should NEVER be used again!!!
setwd("D:/Arabidopsis Arrays")
#find trans markers with plots
mbp <- read.table("refined map/map.physical.c15o15f.txt", row.names=1, header=T)
mcm <- read.table("refined map/map.txt", row.names=1, header=T)
totlength <- c(0, 34964571, 22037565, 25499034, 20862711, 31270811)

addlength_cM <- 0
for(chr in 1:5){
  mbp[which(mbp[,1]==chr),2] <- mbp[which(mbp[,1]==chr),2] + sum(totlength[1:chr])
  
  addlength_cM <- c(addlength_cM, addlength_cM[chr] + max(mcm[which(mcm[,1]==chr),2]))
  mcm[which(mcm[,1]==chr),2] <- mcm[which(mcm[,1]==chr),2] + addlength_cM[chr] + 25*(chr-1)
}
#addlength_cM(0, 123.2792, 196.7278, 289.9451, 391.9624, 505.9731)


#show all chrs within one plot
png(filename = "chrall_bp-cM.png", width = 1024, height = 1024, bg = "white")
plot(mcm[,2], mbp[,2], col=mbp[,1], pch=20, main="Before Remove Markers", xlab="cM", ylab="bp")
dev.off()


#replace trans markers' bp position with NA, and make a new map file
mbp[which(mbp[,1] != mcm[,1]),2] <- NA
mbp["27610",2] <- NA
mbp["126099",2] <- NA
mbp["329114",2] <- NA
write.table(mbp, file="refined map/map_bp_remove_nogap.txt", sep="\t")
mbp_new <- read.table("refined map/map_bp_remove_nogap.txt", row.names=1, header=T)

#show all chrs within one plot
png(filename = "chrall_bp-cM_remove.png", width = 1024, height = 1024, bg = "white")
plot(mcm[,2], mbp_new[,2], col=mbp_new[,1], pch=20, main="After Remove Markers", xlab="cM", ylab="bp")
dev.off()
