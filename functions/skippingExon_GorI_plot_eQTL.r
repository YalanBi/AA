#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 26-07-2013
# first written: 14-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#************************************************************* plot exp part *************************************************************#
#plot 4 env separately in 4 panel
setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]
#load genotype file
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)
#load map file---markers' Morgan position
mpos <- read.table("refined map/map.txt")


#direction selection
probesDir <- function(exp_data=rawexp){
  if(unique(exp_data[ ,"strand"]) == "sense"){
    direction_id <- which(exp_data[ ,"direction"] == "reverse")
  }
  if(unique(exp_data[ ,"strand"]) == "complement"){
    direction_id <- which(exp_data[ ,"direction"] == "forward")
  }
  return(direction_id)
}


dffCol <- function(border=FALSE, x){
  if(!border){
    if(x >= 0) return("red")
    else return("green")
  } else{
    if(x >= 0) return("red4")
    else return("green4")
  }
}


peakDetect <- function(data=testQTL[p, ], cutoff, verbose=FALSE){
  if(verbose) cat("cutoff =", cutoff, "\n")
  
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

########## for marker, use Morgan ##########
totlength <- 0
lengthsx <- -12.5
msumlength <- NULL
for(x in 1:5){
  #msumlength <- each marker's Morgan position
  msumlength <- c(msumlength, mpos[which(mpos[ ,1] == x), 2]+totlength+25*(x-1))
  #totlength <- total length until current chr(Morgan)
  totlength <- totlength+max(mpos[which(mpos[ ,1] == x), 2])
  #lengthsx <- c(the max position of each chr separately(Morgan))
  lengthsx <- c(lengthsx, totlength+25*x-12.5)
}

########## for gene, use bp ##########
chr1 <- 34964571
chr2 <- 22037565
chr3 <- 25499034
chr4 <- 20862711
chr5 <- 31270811
#lengthsy <- c(the max position of each chr separately(bp))
lengthsy <- c(chr1, chr2, chr3, chr4, chr5)

ploteQTL <- function(chr, ind_tu, exonID, rawexp, cutoff=8, verbose=FALSE){
  qtl <- read.table(paste0("Data/FullModel/chr", chr, "_norm_hf_cor_FM/", filename, "_FM_QTL.txt"), row.names=1, header=TRUE)

  plot(c(0.5, nrow(qtl)+0.5), c(0, 610), xaxt='n', xlab="", yaxt='n', ylab="eQTL", cex.axis=1, cex.lab=1.5, las=1, mgp=c(3,1,0), t="n", bg="white")
  for(p in 1:nrow(qtl)){
    if(!p %in% ind_tu) rect((p-0.5), -66, (p+0.5), 666, col=grey(0.85), border = "transparent")#intron probes <- grey background
  }
  
  abline(h=lengthsx[2:5], col="burlywood3", lty=4, lwd=3)
  #axis(2, at=1:5, labels=paste("chr", c("I", "II", "III", "IV", "V"),sep=""), cex.axis=1, las=2, tck=-0.035, mgp=c(2,0.5,0))
  genepos <- median(rawexp[ ,"bp"])/lengthsy[chr]*(lengthsx[chr+1]-lengthsx[1]-25)+lengthsx[chr]+12.5
  abline(h=genepos, col="burlywood3", lty=3, lwd=2)
  #points(-0.5, genepos, pch=">", cex=3, col="burlywood3")
  
  for(p in 1:nrow(qtl)){
    if(p %in% exonID){
      peakList <- peakDetect(data=unlist(qtl[p, ]), cutoff=cutoff)
      
      for(m in 1:716){
        if(peakList[m] == 1){
          rect((p-0.2), (msumlength[m]-1.75), (p+0.2), (msumlength[m]+1.75), col="orange", border="transparent")
          #points(p, msumlength[m], pch=15, cex=1, col="orange")
        }
      }
      for(m in 1:716){
        if(peakList[m] == 2){
          rect((p-0.5), (msumlength[m]-0.5), (p+0.5), (msumlength[m]+0.5), col="chocolate4", border="transparent")
        }
      }
    }
  }
  box();
}


plotExpEnvSep <- function(filename, toTest, testThre, fnPart, m){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  newexp <- rawexp[ ,17:164]
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  ind_tu <- grep("tu", rawexp[ ,"tu"])
  
  dffList <- vector("list", 4)
  for(env in 1:4){
    ind_env <- which(as.numeric(menvironment) == env)
    envGroup1 <- ind_env[ind_env %in% which(geno[,m] == 1)]
    envGroup2 <- ind_env[ind_env %in% which(geno[,m] == 2)]
    if(length(envGroup1) > 0 && length(envGroup2) > 0){
      dff <- NULL
      for(p in probes_dir){
        dff <- c(dff, median(unlist(newexp[p, envGroup1]))-median(unlist(newexp[p, envGroup2])))
      }
      dffList[[env]] <- dff
    }
  }
  
  
  for(env in 1:4){
    ind_env <- which(as.numeric(menvironment) == env)
    envGroup1 <- ind_env[ind_env %in% which(geno[,m] == 1)]
    envGroup2 <- ind_env[ind_env %in% which(geno[,m] == 2)]
    
    if(env == 1){
      par(fig=c(0, 1, 67/64-0.1875*env, 19/16-0.1875*env), oma=c(5, 3, 5, 0.5), mar=c(0, 4, 0, 2))
    }else par(fig=c(0, 1, 67/64-0.1875*env, 19/16-0.1875*env), oma=c(5, 3, 5, 0.5), mar=c(0, 4, 0, 2), new=TRUE)
    plot(c(0.5, nrow(rawexp)+0.5), c(min(rawexp[ ,17:164])-0.2, max(rawexp[ ,17:164])+0.2), xaxt='n', xlab="", ylab=levels(menvironment)[env], cex.axis=1, cex.lab=1.5, col.lab=env+1, las=1, mgp=c(2.25, 0.5, 0), tck=-0.017, t="n")
    if(env == 1){
      title(main=paste0(filename, "_", toTest), sub=paste0("@marker", m), cex.main=2.5, cex.sub=1.25, xlab="", mgp=c(3, 0.5, 0), cex.lab=1.5, outer=TRUE)
      title(ylab="Expression Intensity", mgp=c(1, 0.5, 0), cex.lab=1.5, outer=TRUE)
    }
    
    for(p in 1:nrow(rawexp)){
      if(!p %in% ind_tu) rect((p-0.5), -3, (p+0.5), 13, col=grey(0.85), border="transparent")#background for introns
      if(p %in% probes_dir){
        points(rep(p-0.05, length(envGroup1)), newexp[p, envGroup1], t='p', col="skyblue3", pch=20, cex=0.75)
        points(rep(p+0.05, length(envGroup2)), newexp[p, envGroup2], t='p', col="hotpink3", pch=20, cex=0.75)
        
        if(length(envGroup1) > 0) points(p, median(unlist(newexp[p, envGroup1])), col="skyblue4", pch='*', cex=2)
        if(length(envGroup2) > 0) points(p, median(unlist(newexp[p, envGroup2])), col="hotpink4", pch='*', cex=2)
      }
    }
    
    for(tu in uniqueExon){
      ind <- which(rawexp[ ,"tu"] == tu)
      if(env == 1) axis(3, at=mean(ind), labels=tu, mgp=c(2.25, 0.5, 0), tick=FALSE, line=0)
      #if(any(exonID %in% ind) && median(unlist(newexp[exonID[exonID %in% ind], ])) >= 5){
        #lines(c(min(ind)-0.5, max(ind)+0.5), rep(median(unlist(newexp[exonID[exonID %in% ind], envGroup1])), 2), col="skyblue3", lty=2, lwd=1.5)
        #lines(c(min(ind)-0.5, max(ind)+0.5), rep(median(unlist(newexp[exonID[exonID %in% ind], envGroup2])), 2), col="hotpink3", lty=2, lwd=1.5)
      #}
      
      if(any(ind %in% fnPart[fnPart[ ,3] == m, 2])){
        gseLastProbe <- ind[ind %in% fnPart[fnPart[ ,3] == m, 2]]
        e <- which(fnPart[ ,2] == gseLastProbe)
        if((fnPart[e, 2*env+2]-testThre)*(fnPart[e, 2*env+3]-testThre) <= 0){
          if(fnPart[e, 2*env+2] >= testThre){
            text(x=mean(ind), y=max(newexp)+0.1, labels="gt1 lower", col="skyblue3", cex=0.8)
            #lines(c(min(ind)-0.5, max(ind)+0.5), rep(median(unlist(newexp[exonID[exonID %in% ind], envGroup1])), 2), col="skyblue3", lty=1, lwd=2)
          } #else lines(c(min(ind)-0.5, max(ind)+0.5), rep(median(unlist(newexp[exonID[exonID %in% ind], envGroup1])), 2), col="skyblue3", lty=4, lwd=2)
          if(fnPart[e, 2*env+3] >= testThre){
            text(x=mean(ind), y=max(newexp), labels="gt2 lower", col="hotpink3", cex=0.8)
            #lines(c(min(ind)-0.5, max(ind)+0.5), rep(median(unlist(newexp[exonID[exonID %in% ind], envGroup2])), 2), col="hotpink3", lty=1, lwd=2)
          } #else lines(c(min(ind)-0.5, max(ind)+0.5), rep(median(unlist(newexp[exonID[exonID %in% ind], envGroup2])), 2), col="hotpink3", lty=4, lwd=2)
          if(fnPart[e, 2*env+2] == -1) text(x=mean(ind), y=max(newexp)-0.25, labels="only gt1", cex=0.8)
          if(fnPart[e, 2*env+3] == -2) text(x=mean(ind), y=max(newexp)-0.25, labels="only gt2", cex=0.8)
        } #else{
          #lines(c(min(ind)-0.5, max(ind)+0.5), rep(median(unlist(newexp[exonID[exonID %in% ind], envGroup1])), 2), col="skyblue3", lty=4, lwd=1.5)
          #lines(c(min(ind)-0.5, max(ind)+0.5), rep(median(unlist(newexp[exonID[exonID %in% ind], envGroup2])), 2), col="hotpink3", lty=4, lwd=1.5)
        #}
      }
    }
    box()
    
    par(fig=c(0, 1, 1-0.1875*env, 67/64-0.1875*env), oma=c(5, 3, 5, 0.5),  mar=c(0, 4, 0, 2), new=TRUE)
    if(!is.null(dffList[[env]])){
      plot(c(0.5, nrow(rawexp)+0.5), c(min(unlist(dffList))-0.05, max(unlist(dffList))+0.05), xaxt='n', xlab="", ylab="", cex.axis=1, cex.lab=1.5, col.lab=env, las=1, mgp=c(2.25, 0.5, 0), tck=-0.017, t="n")
      cnt <- 1
      for(p in 1:nrow(rawexp)){
        if(!p %in% ind_tu) rect((p-0.5), -13, (p+0.5), 13, col=grey(0.85), border="transparent")#background for introns
        if(p %in% probes_dir){
          rect(p-0.5, 0, p+0.5, dffList[[env]][cnt], col=dffCol(border=FALSE, dffList[[env]][cnt]), border=dffCol(border=TRUE, dffList[[env]][cnt]))
          cnt <- cnt + 1
        }
      }
    }
    box()
  }
  par(fig=c(0, 1, 0, 0.25), oma=c(5, 3, 5, 0.5),  mar=c(0, 4, 0, 2), new=TRUE)
  ploteQTL(chr, ind_tu, exonID, rawexp, cutoff=8)
  axis(1, at=probes_dir, labels=row.names(rawexp)[probes_dir], mgp=c(2.25, 0.5, 0), cex.axis=1, las=2, tck=0.02)
}
#plotExpEnvSep("AT1G01800", toTest="QTL", testThre=4.50, fnPart, m=2)


#plot all genes having GSE
testThre=4.50
for(chr in 1:5){
  load(file=paste0("Data/geneticsAS/GSE_chr", chr, "_wt_p3.Rdata"))
  genenames <- gseGeneList[[5]]
  for(filename in genenames){
    cat(filename, "starts plotting ...\n")
    
    gsechr <- read.table(paste0("Data/geneticsAS/GSE_chr", chr, "_wt_p3.txt"), row.names=NULL)
    fnPart <- gsechr[grepl(filename, gsechr[ ,1]), ]
    
    sigMarker <- NULL
    for(e in 1:nrow(fnPart)){
      res <- NULL 
      for(env in 1:4){
        if((fnPart[e, 2*env+2]-testThre)*(fnPart[e, 2*env+3]-testThre) <= 0) res <- c(res, 1)
        else res <- c(res, 0)
      }
      if(sum(res) > 0) sigMarker <- c(sigMarker, fnPart[e, 3])
    }
    markerRange <- unique(sigMarker)
    
    for(m in markerRange){
      png(filename = paste0("Data/geneticsAS/Plot/GSE/", filename, "_GSE_m", m, "_wt.png"), width = 960, height = 1728, bg = "white")
      plotExpEnvSep(filename, toTest="QTL", testThre=4.50, fnPart, m)
      dev.off()
    }
  }
}
