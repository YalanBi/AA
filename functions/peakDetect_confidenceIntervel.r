#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 15-04-2013
# first written: 15-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")


#*********************************************************** peak detection part ***********************************************************#
peakDetect <- function(data = qtl[p, ], cutoff = 4){
  #cat("cutoff =", cutoff, "\n")
  
  chrEdges <- c(1, 162, 263, 403, 527)
  
  above   <- as.numeric(data > cutoff)
  below   <- -as.numeric((data < (-cutoff)))
  pattern <- above+below

  inPeek  <- 0  # Are we in a peak ?
  top     <- 0  # What was our previous TOP score ?
  for(x in 1:length(pattern)){
    if(x %in% chrEdges){
      if(inPeek == -1) pattern[topx] <- -2
      if(inPeek == 1) pattern[topx] <- 2
      inPeek <- 0
    }
    #cat(x, pattern[x], data[x], inPeek,"\n")
    if(pattern[x] == 1){
      if(inPeek == 0){
        inPeek <- 1
        top  <- data[x]
        topx <- x
      }else if(inPeek == 1){ # We're in a peak already check this one is the 'top'
        if(top < data[x]){ top <- data[x] ; topx <- x }
        if(top > data[x]){ pattern[topx] <- 2; top <- data[x]; } # Going down, might be a double peek
      }else{ stop(paste("Positive peak to negative peak at", x)) }
    }else if(pattern[x] == -1){
      if(inPeek == 0){
        inPeek <- (-1)
        top  <- data[x]
        topx <- x
      }else if(inPeek == -1){ # We're in a peak already check this one is the 'top'
        if(top > data[x]){ top <- data[x]; topx <- x } # Its a peak
        if(top < data[x]){ pattern[topx] <- -2; top <- data[x]; } # Going down, might be a double peek
      }else{ stop(paste("Negative peak to positive peak at", x)) }
    }else if(pattern[x] == 0){
      if(inPeek == 1){ inPeek = 0; pattern[topx] <- 2 }
      if(inPeek == -1){ inPeek = 0; pattern[topx] <- -2 }
    }
  }
  pattern
}


#*************************************************************** plot part ***************************************************************#
################################################################# panel 1 #################################################################
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]

#plot for expression intensity
makePlot_Exp <- function(filename, rawexp, newexp, nprobes, ind_tu){
  plot(c(0.5, nprobes+0.5), c(round(min(newexp)-0.5), round(max(newexp)+0.5)), xaxt='n', xlab="", ylab="Intensity", cex.axis=1, cex.lab=1.5, cex.main=2, las=1, mgp=c(3,1,0), tck=-0.017, t="n", main=gsub(".txt", "", filename), bg="white")
  s <- 1
  for(p in 1:nprobes){
    #background for introns
    if(!p %in% ind_tu){
      rect((p-0.5), -3, (p+0.5), max(newexp)* 2.1, col=grey(0.85), border = "transparent")
    }
    
    points(rep(p,148) + 0.06 * as.numeric(menvironment) - 0.15, newexp[p, ], t='p', col=menvironment, pch=20, cex=0.7)
    if(rawexp[p,"tu"] != rawexp[s,"tu"]){
      lines(c(s-0.5, p-0.5), c(mean(as.matrix(newexp[s:(p-1),])), mean(as.matrix(newexp[s:(p-1),]))))
      s <- p
    }
  }
  lines(c(s-0.5, p+0.5), c(mean(as.matrix(newexp[s:p,])), mean(as.matrix(newexp[s:p,]))))
  box();
}


################################################################# panel 2 #################################################################
#load map file---markers' Morgan position
mpos <- read.table("refined map/map.txt")

########## for marker, use Morgan ##########
totlength <- 0
lengthsx <- -12.5
msumlength <- NULL
for(x in 1:5){
  #msumlength <- each marker's Morgan position
  msumlength <- c(msumlength, mpos[which(mpos[,1]==x),2] + totlength + 25*(x-1))
  #totlength <- total length until current chr(Morgan)
  totlength = totlength + max(mpos[which(mpos[,1]==x),2])
  #lengthsx <- c(the max position of each chr separately(Morgan))
  lengthsx <- c(lengthsx, totlength + 25*x - 12.5)
}

########## for gene, use bp ##########
chr1 <- 34964571
chr2 <- 22037565
chr3 <- 25499034
chr4 <- 20862711
chr5 <- 31270811
#lengthsy <- c(the max position of each chr separately(bp))
lengthsy <- c(chr1, chr2, chr3, chr4, chr5)

makePlot_eQTL <- function(chr, nprobes, ind_tu, rawexp, qtl, cutoff){
  plot(c(0.5, nprobes+0.5),c(0, 610), xaxt='n', xlab="", yaxt='n', ylab="eQTL", cex.axis=1, cex.lab=1.5, las=1, mgp=c(3,1,0), t="n", bg="white")
  
  for(p in 1:nprobes){
    if(!p %in% ind_tu){
      rect((p-0.5), -66, (p+0.5), 666, col=grey(0.85), border = "transparent")
    }
  }
  
  #cat("cutoff =", cutoff, "\n")
  
  genepos <- median(rawexp[,"bp"])/lengthsy[chr]*(lengthsx[chr+1]-lengthsx[1]-25)+lengthsx[chr]+12.5
  abline(h = genepos, col="burlywood3", lty=3, lwd=2)
  points(-0.5, genepos, pch=">", cex=3, col="burlywood3")
  
    for(p in 1:nprobes){
    peakList <- peakDetect(data = unlist(qtl[p, ]), cutoff = cutoff)
    
    for(m in 1:716){
      if(peakList[m] == 1){
        rect((p-0.2),(msumlength[m]-1.75),(p+0.2),(msumlength[m]+1.75), col="orange", border="transparent")
        #points(p, msumlength[m], pch=15, cex=1, col="orange")
      }
    }
    for(m in 1:716){
      if(peakList[m] == 2){
        rect((p-0.5),(msumlength[m]-0.5),(p+0.5),(msumlength[m]+0.5), col="chocolate4", border="transparent")
      }
    }
  }
  abline(h = lengthsx[2:5], col="burlywood3", lty=4, lwd=3)
  #axis(2, at=1:5, labels=paste("chr", c("I", "II", "III", "IV", "V"),sep=""), cex.axis=1, las=2, tck=-0.035, mgp=c(2,0.5,0))
  
  box();
}


################################################################# panel 3 #################################################################
#means of individuals under each environment for every probe
envTtest2 <- function(newexp, menvironment, p){
  env1 <- mean(as.numeric(newexp[p, which(as.numeric(menvironment)==1)]))
  env2 <- mean(as.numeric(newexp[p, which(as.numeric(menvironment)==2)]))
  env3 <- mean(as.numeric(newexp[p, which(as.numeric(menvironment)==3)]))
  env4 <- mean(as.numeric(newexp[p, which(as.numeric(menvironment)==4)]))
  #lgp <- -log10(c(t.test(env1, mu=mean(c(env2,env3,env4)))$p.value,t.test(env2, mu=mean(c(env1,env3,env4)))$p.value,t.test(env3, mu=mean(c(env1,env2,env4)))$p.value,t.test(env4, mu=mean(c(env1,env2,env3)))$p.value))
  #env_mean <- mean(c(env1,env2,env3,env4))
  #lgp <- -log10(c(t.test(env1, mu=env_mean)$p.value,t.test(env2, mu=env_mean)$p.value,t.test(env3, mu=env_mean)$p.value,t.test(env4, mu=env_mean)$p.value))
  return(c(env1,env2,env3,env4))
}

#load color package and design color pattern
library("colorspace", lib.loc="C:/R/win-library/2.15")
cols <- diverge_hcl(31, h=c(195,330), c = 95, l = c(20, 90), power = 1.25)

#plot of environment comparison
makePlot_Env <- function(nprobes, newexp){
  plot(c(0.5, nprobes+0.5),c(0.5,4.5), xaxt='n', xlab="", yaxt='n', ylab="Env", cex.axis=1, cex.lab=1.5, las=1, mgp=c(3,1,0), t="n", bg="white")
  axis(1, at=1:nprobes, labels=row.names(newexp), cex.axis=1, las=2, tck=0.035, mgp=c(2,1,0))
  axis(2, at=1, labels="6H", cex.axis=1, col.axis='black', las=2, tck=-0.02, mgp=c(2,0.5,0))
  axis(2, at=2, labels="DA", cex.axis=1, col.axis='red', las=2, tck=-0.02, mgp=c(2,0.5,0))
  axis(2, at=3, labels="DF", cex.axis=1, col.axis='green', las=2, tck=-0.02, mgp=c(2,0.5,0))
  axis(2, at=4, labels="RP", cex.axis=1, col.axis='blue', las=2, tck=-0.02, mgp=c(2,0.5,0))
  for(p in 1:nprobes){
    #pForCol <- round(envTtest(newexp,env,p)-0.5)+1
    pForCol <- ((envTtest2(newexp, menvironment, p) - mean(as.numeric(newexp[p,])))*7.5) + 16
    pForCol[pForCol > 31] <- 31
    pForCol[pForCol < 1] <- 1
    #cat(pForCol,"\n")
    rect((p-0.5),0.6,(p+0.5),1.4,col=cols[pForCol[1]],border = "transparent")
    rect((p-0.5),1.6,(p+0.5),2.4,col=cols[pForCol[2]],border = "transparent")
    rect((p-0.5),2.6,(p+0.5),3.4,col=cols[pForCol[3]],border = "transparent")
    rect((p-0.5),3.6,(p+0.5),4.4,col=cols[pForCol[4]],border = "transparent")
  }
  box();
}

################################################################# make 3 in 1 #################################################################
singlePlot <- function(chr, filename, cutoff = 4){
  rawexp <- read.table(paste0(location, filename), row.names=1, header=T)
  qtl <- read.table(paste0(location_FM, gsub(".txt", "_FM_QTL.txt", filename)), row.names=1, header=T)
  ind_tu <- grep("tu", rawexp[,"tu"])
  newexp <- rawexp[,17:164]
  nprobes <- nrow(qtl)
  
  st <- proc.time()
  
  par(fig=c(0,1,0.7,1), mar=c(0,5,5,2), oma=c(0,0,0,0.5), bg="white")
  makePlot_Exp(filename, rawexp, newexp, nprobes, ind_tu)
  par(fig=c(0,1,0.28,0.7), mar=c(0,5,0,2), oma=c(0,0,0,0.5), new=T)
  makePlot_eQTL(chr, nprobes, ind_tu, rawexp, qtl, cutoff)
  par(fig=c(0,1,0,0.28), mar=c(5,5,0,2), oma=c(0,0,0,0.5), new=T)
  makePlot_Env(nprobes, newexp)
  
  et <- proc.time()
  cat("Done with", filename, "after:",(et-st)[3],"secs\n")
}

for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  location_FM <- paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/")
  
  genenames <- gsub("_FM_QTL.txt", ".txt", dir(location_FM)[which(grepl("QTL", dir(location_FM)))])
  
  #filename="AT1G01010"
  for(filename in genenames){
    tiff(file = paste(location_FM, gsub(".txt", "_CI.tiff", filename), sep=""), bg="white", width=1536, height=1536, pointsize = 18, compression="none", res=48)
    singlePlot(chr, filename, cutoff = 4)
    dev.off()
  }
}
