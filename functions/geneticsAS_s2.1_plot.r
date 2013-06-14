#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 14-06-2013
# first written: 14-06-2013
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
dffCol <- function(border = FALSE, x){
  if(!border){
    if(x >= 0) return("red")
    else return("green")
  } else{
    if(x >= 0) return("red4")
    else return("green4")
  }
}


peakDetect <- function(data = testQTL[p, ], cutoff = 4){
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

makePlot_eQTL <- function(chr, nprobes, ind_tu, exonID, rawexp, qtl = testQTL, cutoff){
  plot(c(0.5, nprobes+0.5), c(0, 610), xaxt='n', xlab="", yaxt='n', ylab="eQTL", cex.axis=1, cex.lab=1.5, las=1, mgp=c(3,1,0), t="n", bg="white")
  
  for(p in 1:nprobes){
    if(!p %in% ind_tu){
      rect((p-0.5), -66, (p+0.5), 666, col=grey(0.85), border = "transparent")
    }
  }
  
  #cat("cutoff =", cutoff, "\n")
  
  genepos <- median(rawexp[,"bp"])/lengthsy[chr]*(lengthsx[chr+1]-lengthsx[1]-25)+lengthsx[chr]+12.5
  abline(h = genepos, col="burlywood3", lty=3, lwd=2)
  #points(-0.5, genepos, pch=">", cex=3, col="burlywood3")
  
  for(p in 1:nprobes){
    if(p %in% exonID){
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
  }
  abline(h = lengthsx[2:5], col="burlywood3", lty=4, lwd=3)
  #axis(2, at=1:5, labels=paste("chr", c("I", "II", "III", "IV", "V"),sep=""), cex.axis=1, las=2, tck=-0.035, mgp=c(2,0.5,0))
  
  box();
}


plotExpByGtEnvSep <- function(rawexp, m, ind_env, use = median, verbose = FALSE){
  newexp <- rawexp[ ,17:164]
  probes_dir <- probesDir(rawexp)
  ind_tu <- grep("tu", rawexp[ ,"tu"])
  
  geno1 <- which(geno[,m] == 1)
  geno2 <- which(geno[,m] == 2)
  envGroup1 <- ind_env[ind_env %in% geno1]
  envGroup2 <- ind_env[ind_env %in% geno2]
  
  if(verbose) cat(use, "plotExp\n")
  for(p in 1:nrow(rawexp)){
    #background for introns
    if(!p %in% ind_tu) rect((p - 0.5), -3, (p + 0.5), 13, col = grey(0.85), border = "transparent")
    if(p %in% probes_dir){
      points(rep(p-0.05, length(envGroup1)), newexp[p,envGroup1], t = 'p', col = "skyblue3", pch = 20, cex = 0.75)
      points(rep(p+0.05, length(envGroup2)), newexp[p,envGroup2], t = 'p', col = "hotpink3", pch = 20, cex = 0.75)
      
      points(p, apply(newexp[p,envGroup1], 1, use), col = "skyblue4", pch = '*', cex = 1.5)
      points(p, apply(newexp[p,envGroup2], 1, use), col = "hotpink4", pch = '*', cex = 1.5)
    }
  }
  box()
}

dffBtwGt <- function(rawexp, m, ind = ind_env, use = median, verbose = FALSE){
  newexp <- rawexp[ ,17:164]
  probes_dir <- probesDir(rawexp)
  
  geno1 <- which(geno[,m] == 1)
  geno2 <- which(geno[,m] == 2)
  envGroup1 <- ind[ind %in% geno1]
  envGroup2 <- ind[ind %in% geno2]
  
  if(verbose) cat(use, "dffBteGt\n")
  dff <- NULL
  if(length(envGroup1) > 0 && length(envGroup2) > 0){
    for(p in 1:nrow(rawexp)){
      if(p %in% probes_dir){
        dff <- c(dff, apply(newexp[p,envGroup1], 1, use) - apply(newexp[p,envGroup2], 1, use))
      }
    }
    return(dff)
  } else return(NULL)
}

plotExpDffBtwGtEnvSep <- function(rawexp, dff = dffList[[env]]){
  probes_dir <- probesDir(rawexp)
  ind_tu <- grep("tu", rawexp[ ,"tu"])
  
  dff_cnt <- 1
  for(p in 1:nrow(rawexp)){
    #background for introns
    if(!p %in% ind_tu) rect((p - 0.5), -3, (p + 0.5), 13, col = grey(0.85), border = "transparent")
    if(p %in% probes_dir){
      rect(p-0.5, 0, p+0.5, dff[dff_cnt], col=dffCol(border = FALSE, dff[dff_cnt]), border=dffCol(border = TRUE, dff[dff_cnt]))
      dff_cnt <- dff_cnt + 1
    }
  }
  box()
}



    for(tu in asTU){
      if(env == 1) axis(3, at = mean(which(rawexp[ ,"tu"] == tu)), labels = tu, mgp=c(2.25, 0.5, 0), tick = FALSE, line = 0)
      if(GASmatrix[grep(paste0(filename, "_", tu), rownames(GASmatrix), value=T), env][1] >= gas_thre) text(x = mean(which(rawexp[ ,"tu"] == tu)), y = max(newexp)+0.1, labels="gt1 lower", col="skyblue3", cex = 0.8)
      if(GASmatrix[grep(paste0(filename, "_", tu), rownames(GASmatrix), value=T), env][2] >= gas_thre) text(x = mean(which(rawexp[ ,"tu"] == tu)), y = max(newexp), labels="gt2 lower", col="hotpink3", cex = 0.8)
      if(GASmatrix[grep(paste0(filename, "_", tu), rownames(GASmatrix), value=T), env][1] == -1) text(x = mean(which(rawexp[ ,"tu"] == tu)), y = max(newexp)-0.25, labels="only gt1", cex = 0.8)
      if(GASmatrix[grep(paste0(filename, "_", tu), rownames(GASmatrix), value=T), env][2] == -2) text(x = mean(which(rawexp[ ,"tu"] == tu)), y = max(newexp)-0.25, labels="only gt2", cex = 0.8)
    }


plotExpEnvSep <- function(filename, toTest, m, ...){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  
  rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  testQTL <- read.table(paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/", filename, "_FM_", toTest, ".txt"), row.names=1, header=T)
  
  probes_dir <- probesDir(rawexp)
  #cat(filename, "\nprobeDir:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #cat("exons of right direction:", exonID, "\n")
  
  ind_tu <- grep("tu", rawexp[ ,"tu"])
  
  dffList <- vector("list", 4)
  for(env in 1:4){
    ind_env <- which(as.numeric(menvironment) == env)
    dffInEnv <- dffBtwGt(rawexp, m, ind = ind_env, ...)
    if(length(dffInEnv) > 0) dffList[[env]] <- dffInEnv
  }
  
  #par(mfrow = c(4, 1), pty = "m", oma = c(5, 3, 5, 0.5),  mar = c(0, 4, 0, 2))
  for(env in 1:4){
    ind_env <- which(as.numeric(menvironment) == env)
    
    if(env == 1){
      par(fig = c(0, 1, 67/64-0.1875*env, 19/16-0.1875*env), oma = c(5, 3, 5, 0.5), mar = c(0, 4, 0, 2))
    }else par(fig = c(0, 1, 67/64-0.1875*env, 19/16-0.1875*env), oma = c(5, 3, 5, 0.5), mar = c(0, 4, 0, 2), new = TRUE)
    plot(c(0.5, nrow(rawexp) + 0.5), c(min(rawexp[ ,17:164]) - 0.2, max(rawexp[ ,17:164]) + 0.2), xaxt = 'n', xlab = "", ylab = levels(menvironment)[env], cex.axis = 1, cex.lab = 1.5, las = 1, mgp = c(2.25, 0.5, 0), tck = -0.017, t = "n")
    if(env == 1){
      title(main = paste0(filename, "_", toTest), sub = paste0("@marker", m), cex.main = 2.5, cex.sub = 1.25, xlab = "", mgp = c(3, 0.5, 0), cex.lab = 1.5, outer = TRUE)
      title(ylab = "Expression Intensity", mgp = c(1, 0.5, 0), cex.lab = 1.5, outer = TRUE)
    }
    plotExpByGtEnvSep(rawexp, m, ind_env, ...)
    
    par(fig = c(0, 1, 1-0.1875*env, 67/64-0.1875*env), oma = c(5, 3, 5, 0.5),  mar = c(0, 4, 0, 2), new = TRUE)
    if(length(dffList[[env]]) > 0){
      plot(c(0.5, nrow(rawexp) + 0.5), c(min(unlist(dffList)) - 0.05, max(unlist(dffList)) + 0.05), xaxt = 'n', xlab = "", ylab = "", cex.axis = 1, cex.lab = 1.5, col.lab = env, las = 1, mgp = c(2.25, 0.5, 0), tck = -0.017, t = "n")
      plotExpDffBtwGtEnvSep(rawexp, dff = dffList[[env]])
    }
    box()
  }
  par(fig = c(0, 1, 0, 0.25), oma = c(5, 3, 5, 0.5),  mar = c(0, 4, 0, 2), new = TRUE)
  makePlot_eQTL(chr, nprobes=nrow(rawexp), ind_tu, exonID, rawexp, qtl = testQTL, cutoff = 8)
  axis(1, at = probes_dir, labels = row.names(rawexp)[probes_dir], mgp=c(2.25, 0.5, 0), cex.axis = 1, las = 2, tck = 0.02)
}
plotExpEnvSep("AT1G01800", toTest="QTL", m=2, use="median", verbose=T)



#ttestExp
toTest <- "QTL"#"Int"
GASmatrix <- NULL; gas_thre = -log10(0.05/(616*4))#-log10(0.05/(39*4))
for(chr in 1:5){
  GASmatrix <- rbind(GASmatrix, read.table(paste0("Data/countQTL/main", toTest, "_chr", chr, "_ttestExp.txt"), row.names=1, header=F))
}
allGenenames <- unique(unlist(lapply(strsplit(rownames(GASmatrix), "_"), "[[", 1)))
potentialGAS <- NULL
for(filename in allGenenames){
  for(tu in unique(unlist(lapply(strsplit(grep(filename, rownames(GASmatrix), value=T), "_"), "[[", 2)))){
    for(env in 1:4){
      if(any(GASmatrix[grep(paste0(filename, "_", tu), rownames(GASmatrix), value=T), env] >= gas_thre) && any(GASmatrix[grep(paste0(filename, "_", tu), rownames(GASmatrix), value=T), env] <= gas_thre)){
        potentialGAS <- c(potentialGAS, filename)
      }
    }
  }
}
plotGenenames <- unique(potentialGAS)
#filename(AT1G01010)
for(filename in plotGenenames){
  png(filename = paste0("Data/countQTL/GASPlot/", filename, "_", toTest, "_4s.png"), width = 960, height = 1728, bg = "white")
  plotExpEnvSep(filename, toTest = "QTL", GASmatrix, gas_thre, use = mean)#plotExpEnvSep(filename, toTest = "Int", GASmatrix, gas_thre, use = mean)
  dev.off()
}
#plot for ttestExp
(filename, toTest = "QTL"){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  GASmatrix <- read.table(paste0("Data/countQTL/main", toTest, "_chr", chr, "_ttestExp.txt"), row.names=1, header=F)
  rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  testQTL <- read.table(paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/", filename, "_FM_", toTest, ".txt"), row.names=1, header=T)
  
  probes_dir <- probesDir(rawexp)
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  
  asTU <- unique(unlist(lapply(strsplit(grep(filename, rownames(GASmatrix), value=TRUE), "_"), "[[", 2)))
  
  mlist <- NULL
  for(tu in asTU){
    tuID <- exonID[rawexp[exonID, "tu"] == tu]
    mlist <- c(mlist, which.max(apply(testQTL[tuID,], 2, sum)))
    cat(tu, "m", which.max(apply(testQTL, 2, sum)), "\n")
  }
  cat("whole gene m", which.max(apply(testQTL, 2, sum)), "\n")
  for(nplot in 1:length(unique(mlist))){
    m <- unique(mlist)[nplot]
    geno1 <- which(geno[,m] == 1)
    geno2 <- which(geno[,m] == 2)
    plotExpEnvSep(filename, toTest="QTL", m, use="median")
  }
}

