#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 16-07-2013
# first written: 16-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#*************************************************** this is the final version for plotting alternative 5'/3' site ^_^ ***************************************************#
#**************************************************************** testing algorithm: Wilcox.test / ANOVA! ****************************************************************#

setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]

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

plot_53AS <- function(filename, goal, whichTest, ){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  newexp <- rawexp[ ,17:164]
  probes_dir <- probesDir(rawexp)
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  ind_tu <- grep("tu", rawexp[ ,"tu"])
  nprobes <- nrow(rawexp)
  
  testAS <- read.table(paste0("Data/AS/", goal, "_chr", chr, "_", whichTest, "_less.txt"), row.names=NULL)
  fnPart <- testAS[testAS[ ,1] %in% filename, ]
  
  par(mfrow=c(4, 1), pty="m", oma=c(5, 3, 5, 0.5), mar=c(0, 4, 0, 2))
  for(env in 1:4){
    ind_env <- which(as.numeric(menvironment) == env)
    plot(c(0.5, nprobes+0.5), c(min(newexp)-0.2, max(newexp)+0.2), xaxt='n', xlab="", ylab=levels(menvironment)[env], cex.axis=1, cex.lab=1.5, col.lab=env+1, las=1, mgp=c(2.25, 0.5, 0), tck=-0.017, t="n")
    
    if(env == 1){
      title(main=filename, cex.main=2.5, sub=paste0(goal, "_", whichTest), cex.sub=2, xlab="Probes", cex.lab=1.5, mgp=c(3, 0.5, 0), outer=TRUE)
      title(ylab="Expression Intensity", mgp=c(1, 0.5, 0), cex.lab=1.5, outer=TRUE)
      for(tu in uniqueExon){
        axis(3, at=mean(which(rawexp[ ,"tu"] == tu)), labels=tu, mgp=c(2.25, 0.5, 0), tick=FALSE, line=0)
      }
    }
    
    if(goal == "5'AS"){#******************************************************************** HERE! ********************************************************************#
      asThre=5.05
      sepPoint <- fnPart[ ,2]
      rect(0.5, -3, sepPoint+0.5, 13, density = 15, angle = 30, col = "azure2", border = "azure3")
      rect(sepPoint+0.5, -3, max(which(rawexp[ ,"tu"] == uniqueExon[1]))+0.5, 13, density = 15, angle = 150, col = "azure2", border = "azure3")
      
      lines(c(0.5, sepPoint+0.5), rep(mean(unlist(newexp[-(sepPoint:max(which(rawexp[, "tu"] == uniqueExon[1]))), ind_env])), 2), lwd=2)
      if(fnPart[ ,env+2] >= asThre){
        text(x=mean(c(sepPoint, max(which(rawexp[, "tu"] == uniqueExon[1])))), y=max(newexp)+0.1, labels=round(fnPart[e, env+2], digits=1), col=env+1, cex=1)
        lines(c(sepPoint+0.5, max(which(rawexp[, "tu"] == uniqueExon[1]))+0.5), rep(mean(unlist(newexp[sepPoint:max(which(rawexp[, "tu"] == uniqueExon[1])), ind_env])), 2), col=env+1, lwd=2)
      } else{
        if(round(fnPart[e, env+2], digits=2) == asThre) text(x=mean(c(sepPoint, max(which(rawexp[, "tu"] == uniqueExon[1])))), y=max(newexp)+0.1, labels=round(fnPart[e, env+2], digits=1), col=1, cex=1)
        lines(c(sepPoint+0.5, max(which(rawexp[, "tu"] == uniqueExon[1]))+0.5), rep(mean(unlist(newexp[sepPoint:max(which(rawexp[, "tu"] == uniqueExon[1])), ind_env])), 2), col=1, lwd=2)
      }
      lines(c(max(which(rawexp[, "tu"] == uniqueExon[1]))+0.5, nprobes+0.5), rep(mean(unlist(newexp[-(sepPoint:max(which(rawexp[, "tu"] == uniqueExon[1]))),ind_env])), 2), lwd=2)
    }
    
    if(goal == "3'AS"){
      asThre=5.09
    }
    
    
    
    
    for(p in 1:nprobes){
      #background for introns
      if(!p %in% ind_tu){
        rect((p-0.5), -3, (p+0.5), max(newexp[ ,ind_env])*2, col=grey(0.85), border="transparent")
      }
      if(p %in% probes_dir){
        points(rep(p, length(ind_env))+runif(length(ind_env), min=-0.05, max=0.05), newexp[p, ind_env], t='p', pch=20, cex=0.75)
      }
    }
    box()
  }
  axis(1, at=probes_dir, labels=row.names(rawexp)[probes_dir], mgp=c(2.25, 0.5, 0), cex.axis=1, las=2, tck=0.02)
}


#plot all 5'3'AS genes
goal="5'AS"# "3'AS"
whichTest="wt"# "ANOVA"
for(chr in 1:5){
  load(file=paste0("Data/AS/", goal, "_chr", chr, "_", whichTest, ".Rdata"))
  genenames <- asGeneList$mixEnv
  
  for(filename in genenames){
    png(filename = paste0("Data/AS/plot/5'3'AS/", filename, "_", goal, ".png"), width = 960, height = 1728, bg = "white")
    plot_53AS(filename, whichTest, asThre=6.06)
    dev.off()
  }
}

plot_53AS <- function(filename, goal, whichTest, ){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  newexp <- rawexp[ ,17:164]
  
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")

  ind_tu <- grep("tu", rawexp[ ,"tu"])
  nprobes <- nrow(rawexp)
  
  
  testAS <- read.table(paste0("Data/AS/", goal, "_chr", chr, "_", whichTest, "_less.txt"), row.names=NULL)
  
  par(mfrow=c(4, 1), pty="m", oma=c(5, 3, 5, 0.5), mar=c(0, 4, 0, 2))
  for(env in 1:4){
    ind_env <- which(as.numeric(menvironment) == env)
    plot(c(0.5, nprobes+0.5), c(min(newexp)-0.2, max(newexp)+0.2), xaxt='n', xlab="", ylab=levels(menvironment)[env], cex.axis=1, cex.lab=1.5, las=1, mgp=c(2.25, 0.5, 0), tck=-0.017, t="n")
    
    if(env == 1){
      title(main=filename, cex.main=2.5, sub="5'3'AS", cex.sub=2, xlab="Probes", cex.lab=1.5, mgp=c(3, 0.5, 0), outer=TRUE)
      title(ylab="Expression Intensity", mgp=c(1, 0.5, 0), cex.lab=1.5, outer=TRUE)
      for(tu in uniqueExon){
        axis(3, at = mean(which(rawexp[ ,"tu"] == tu)), labels = tu, mgp=c(2.25, 0.5, 0), tick = FALSE, line = 0)
      }
    }
    
    
    for(p in 1:nprobes){
      #background for introns
      if(!p %in% ind_tu){
        rect((p - 0.5), -3, (p + 0.5), max(newexp[ ,ind_env])* 2, col = grey(0.85), border = "transparent")
      }
      if(p %in% probes_dir){
        points(rep(p, length(ind_env)) + runif(length(ind_env), min = -0.05, max = 0.05), newexp[p,ind_env], t = 'p', pch = 20, cex = 0.75)
      }
    }
    box()
  }

}

