#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 10-06-2013
# first written: 10-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#************************************************************* plot exp part *************************************************************#
#plot 4 env separately in 4 panel
setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]
#load genotype file
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)
#load map file---markers' Morgan position
mpos <- read.table("refined map/map.txt")


#direction selection
probesDir <- function(rawexp, minExpression = -1){
  if(unique(rawexp[,"strand"]) == "sense"){
    direction_id <- which(rawexp[, "direction"] == "reverse")
  }
  if(unique(rawexp[,"strand"]) == "complement"){
    direction_id <- which(rawexp[, "direction"] == "forward")
  }
  
  if(minExpression != -1){
    expressed <- which(apply(rawexp[ ,17:164], 1, median) >= minExpression)
    return(direction_id[direction_id %in% expressed])
  }else return(direction_id)
}


genoCol <- function(x){
  x[which(x == 1)] <- "skyblue3"
  x[which(x == 2)] <- "hotpink3"
  return(x)
}
dffCol <- function(x){
  if(x >= 0) return("red")
  else return("green")
}


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

makePlot_eQTL <- function(chr, nprobes, ind_tu, rawexp, qtl = testQTL, cutoff){
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


plotExpEnvSep <- function(filename, toTest = "QTL", use = mean){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  
  rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  newexp <- rawexp[ ,17:164]
  testQTL <- read.table(paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/", filename, "_FM_", toTest, ".txt"), row.names=1, header=T)
  
  probes_dir <- probesDir(rawexp)
  #cat(filename, "\nprobeDir:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #cat("exons of right direction:", exonID, "\n")
  
  #uniqueExon <- all tu names of exon probes
  uniqueExon <- unique(rawexp[exonID,"tu"])
  #cat("tu names:", as.character(uniqueExon), "\n")
  
  ind_tu <- grep("tu", rawexp[ ,"tu"])
  nprobes <- nrow(rawexp)
  
  m <- which.max(apply(testQTL, 2, sum))

  dffList <- vector("list", 4)
  
  #par(mfrow = c(4, 1), pty = "m", oma = c(5, 3, 5, 0.5),  mar = c(0, 4, 0, 2))
  for(env in 1:4){
    ind_env <- which(as.numeric(menvironment) == env)
    
    if(env == 1){
      par(fig = c(0, 1, 17/16-0.1875*env, 19/16-0.1875*env), oma = c(5, 3, 5, 0.5),  mar = c(0, 4, 0, 2))
    }else par(fig = c(0, 1, 17/16-0.1875*env, 19/16-0.1875*env), oma = c(5, 3, 5, 0.5),  mar = c(0, 4, 0, 2), new = TRUE)
    plot(c(0.5, nprobes + 0.5), c(min(newexp) - 0.2, max(newexp) + 0.2), xaxt = 'n', xlab = "", ylab = levels(menvironment)[env], cex.axis = 1, cex.lab = 1.5, las = 1, mgp = c(2.25, 0.5, 0), tck = -0.017, t = "n")
    if(env == 1){
      title(main = filename, cex.main = 2.5, xlab = "", mgp = c(3, 0.5, 0), cex.lab = 1.5, outer = TRUE)
      title(ylab = "Expression Intensity", mgp = c(1, 0.5, 0), cex.lab = 1.5, outer = TRUE)
      axis(3, at = mean(exonID[rawexp[exonID, "tu"] == asTU]), labels = asTU, mgp=c(2.25, 0.5, 0), tick = FALSE, line = 0)
    }
    
    for(p in 1:nprobes){
      #background for introns
      if(!p %in% ind_tu) rect((p - 0.5), -3, (p + 0.5), max(newexp[ ,ind_env])* 2, col = grey(0.85), border = "transparent")
      if(p %in% probes_dir){
        points(rep(p, length(ind_env)) + runif(length(ind_env), min = -0.05, max = 0.05), newexp[p,ind_env], t = 'p', col = genoCol(geno[,m]), pch = 20, cex = 0.75)
        dffList[[env]] <- c(dffList[[env]], apply(newexp[p,ind_env[ind_env %in% which(geno[,m] == 1)]], 1, use)-apply(newexp[p,ind_env[ind_env %in% which(geno[,m] == 2)]], 1, use))
      }
    }
    box()
    
    
    par(fig = c(0, 1, 1-0.1875*env, 17/16-0.1875*env), oma = c(5, 3, 5, 0.5),  mar = c(0, 4, 0, 2), new = TRUE)
    plot(c(0.5, nprobes + 0.5), c(min(unlist(dffList)) - 0.05, max(unlist(dffList)) + 0.05), xaxt = 'n', xlab = "", ylab = "", cex.axis = 1, cex.lab = 1.5, col.lab = env, las = 1, mgp = c(2.25, 0.5, 0), tck = -0.017, t = "n")
    
    dff_cnt <- 1
    for(p in 1:nprobes){
      #background for introns
      if(!p %in% ind_tu) rect((p - 0.5), -3, (p + 0.5), max(newexp[ ,ind_env])* 2, col = grey(0.85), border = "transparent")
      if(p %in% probes_dir){
        rect(p-0.5, 0, p+0.5, dffList[[env]][dff_cnt], col=dffCol(dffList[[env]][dff_cnt]), border="transparent")
        dff_cnt <- dff_cnt + 1
      }
    }
    
    box()
  }
  par(fig = c(0, 1, 0, 0.25), oma = c(5, 3, 5, 0.5),  mar = c(0, 4, 0, 2), new = TRUE)
  makePlot_eQTL(chr, nprobes, ind_tu, rawexp, qtl = testQTL, cutoff = 8)
  axis(1, at = probes_dir, labels = row.names(rawexp)[probes_dir], mgp=c(2.25, 0.5, 0), cex.axis = 1, las = 2, tck = 0.02)
}


GASmatrix <- NULL; gas_thre = -log10(0.05/644)
for(chr in 1:5){
  GASmatrix <- rbind(GASmatrix, read.table(paste0("Data/countQTL/mainQTL_chr", chr, "_ttest.txt"), row.names=1, header=F))
  
  sigTU <- rownames(GASmatrix)[GASmatrix >= gas_thre]
  plotGenenames <- unlist(lapply(strsplit(sigTU, "_"), "[[", 1))
  #filename(AT1G01010)
  for(filename in plotGenenames){
    asTU <- unlist(strsplit(grep(filename, sigTU, value=TRUE), "_"))[2]
    plotcExonExp(chr, filename, ce_threshold = 5.86)
  }
  png(filename = paste0("Data/geneticsAS/plotGeneticsAS/", filename, "_QTL", threshold_qtl, "_Int", threshold_int, "_np", cutoffnProbe, "_4s.png"), width = 960, height = 1728, bg = "white")
  plotExpEnvSep(filename, toTest = "QTL", use = mean)
  dev.off()
}
