#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 26-07-2013
# first written: 26-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#************************************************************* plot exp part *************************************************************#
#plot 4 env separately in 4 panel
setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]
#load genotype file
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)

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
      par(fig=c(0, 1, 17/16-0.25*env, 1.25-0.25*env), oma=c(5, 3, 5, 0.5), mar=c(0, 4, 0, 2))
    }else par(fig=c(0, 1, 17/16-0.25*env, 1.25-0.25*env), oma=c(5, 3, 5, 0.5), mar=c(0, 4, 0, 2), new=TRUE)
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
        if(fnPart[e, 2*env+2] >= testThre && fnPart[e, 2*env+3] >= testThre) text(x=rep(mean(ind), 2), y=c(max(newexp)+0.1, max(newexp)), labels=paste0("both/", round(fnPart[e, (2*env+2):(2*env+3)], digits=1)), col=c("skyblue3", "hotpink3"), cex=0.8)
        else if((fnPart[e, 2*env+2]-testThre)*(fnPart[e, 2*env+3]-testThre) <= 0){
          if(fnPart[e, 2*env+2] >= testThre){
            text(x=mean(ind), y=max(newexp)+0.1, labels=paste0("gt1/", round(fnPart[e, 2*env+2], digits=1)), col="skyblue3", cex=0.8)
            #lines(c(min(ind)-0.5, max(ind)+0.5), rep(median(unlist(newexp[exonID[exonID %in% ind], envGroup1])), 2), col="skyblue3", lty=1, lwd=2)
          } #else lines(c(min(ind)-0.5, max(ind)+0.5), rep(median(unlist(newexp[exonID[exonID %in% ind], envGroup1])), 2), col="skyblue3", lty=4, lwd=2)
          if(fnPart[e, 2*env+3] >= testThre){
            text(x=mean(ind), y=max(newexp), labels=paste0("gt2/", round(fnPart[e, 2*env+3], digits=1)), col="hotpink3", cex=0.8)
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
    
    par(fig=c(0, 1, 1-0.25*env, 17/16-0.25*env), oma=c(5, 3, 5, 0.5),  mar=c(0, 4, 0, 2), new=TRUE)
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


#plot most sig SE
genenames <- NULL
gseMatrix <- NULL
for(chr in 1:5){
  load(file=paste0("Data/geneticsAS/GSE_chr", chr, "_wt_p3.Rdata"))
  genenames <- c(genenames, gseGeneList[[5]])
  if(file.exists(paste0("Data/geneticsAS/GSE_chr", chr, "_wt_p3.txt"))){
    gseMatrix <- rbind(gseMatrix, read.table(paste0("Data/geneticsAS/GSE_chr", chr, "_wt_p3.txt"), row.names=NULL))
  } else cat("chr", chr, "NO test!\n")
}
resList <- NULL
for(e in 1:nrow(gseMatrix)){
  res <- 1
  for(env in 1:4){
    res <- res*abs(gseMatrix[e, 2*env+2]-gseMatrix[e, 2*env+3])
  }
  resList <- c(resList, res)
}
sigOrder <- unique(gseMatrix[order(resList, decreasing=TRUE), 1])
mostSig <- sigOrder[sigOrder %in% genenames]
for(filename in mostSig[1:5]){
  cat(filename, "starts plotting ...\n")
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
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
  cat(markerRange, "\n")
  for(m in markerRange){
    png(filename = paste0("Data/geneticsAS/Plot/GSE/", filename, "_GSE_m", m, "_wt_mS.png"), width = 960, height = 1728, bg = "white")
    plotExpEnvSep(filename, toTest="QTL", testThre=4.50, fnPart, m)
    dev.off()
  }
}
#plot around thre SE
testThre=4.50
genenames <- NULL
gseMatrix <- NULL
for(chr in 1:5){
  load(file=paste0("Data/geneticsAS/GSE_chr", chr, "_wt_p3.Rdata"))
  genenames <- c(genenames, gseGeneList[[5]])
  if(file.exists(paste0("Data/geneticsAS/GSE_chr", chr, "_wt_p3.txt"))){
    gseMatrix <- rbind(gseMatrix, read.table(paste0("Data/geneticsAS/GSE_chr", chr, "_wt_p3.txt"), row.names=NULL))
  } else cat("chr", chr, "NO test!\n")
}
res <- NULL
for(e in 1:nrow(gseMatrix)){
  for(gt in 4:11){
    if(gseMatrix[e, gt] > (testThre-0.1) && gseMatrix[e, gt] < (testThre+0.1)) res <- c(res, gseMatrix[e, 1])
  }
}
aroundThre <- unique(res)[unique(res) %in% genenames]
for(filename in aroundThre[c(1:3, 9)]){
  cat(filename, "starts plotting ...\n")
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
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
  cat(markerRange, "\n")
  for(m in markerRange){
    png(filename = paste0("Data/geneticsAS/Plot/GSE/", filename, "_GSE_m", m, "_wt_aT.png"), width = 960, height = 1728, bg = "white")
    plotExpEnvSep(filename, toTest="QTL", testThre=4.50, fnPart, m)
    dev.off()
  }
}
