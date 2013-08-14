#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 29-07-2013
# first written: 21-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#********************************************** This is used to plot alternative splicing at 5' **********************************************#
#plot 4 env separately in 4 panel
setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]


#direction selection
probesDir <- function(exp_data = rawexp){
  if(unique(exp_data[,"strand"]) == "sense"){
    direction_id <- which(exp_data[, "direction"] == "reverse")
  }
  if(unique(exp_data[,"strand"]) == "complement"){
    direction_id <- which(exp_data[, "direction"] == "forward")
  }
  return(direction_id)
}

plot_3AS <- function(filename, asThre=4.56){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  newexp <- rawexp[ ,17:164]
  probes_dir <- probesDir(rawexp)
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  ind_tu <- grep("tu", rawexp[ ,"tu"])
  
  testAS <- read.table(paste0("Data/AS/3'AS_chr", chr, "_wt_less_p2.txt"), row.names=NULL)
  fnPart <- testAS[testAS[ ,1] %in% filename, ]
  
  par(mfrow=c(4, 1), pty="m", oma=c(5, 3, 5, 0.5), mar=c(0, 4, 0, 2))
  for(env in 1:4){
    ind_env <- which(as.numeric(menvironment) == env)
    plot(c(0.5, nrow(rawexp)+0.5), c(min(newexp)-0.2, max(newexp)+0.2), xaxt='n', xlab="", ylab=levels(menvironment)[env], cex.axis=1, cex.lab=1.5, col.lab=env+1, las=1, mgp=c(2.25, 0.5, 0), tck=-0.017, t="n")
    
    if(env == 1){
      title(main=filename, cex.main=2.5, sub="3'AS by Wilcoxon Test", cex.sub=2, xlab="Probes", cex.lab=1.5, mgp=c(3, 0.5, 0), outer=TRUE)
      title(ylab="Expression Intensity", mgp=c(1, 0.5, 0), cex.lab=1.5, outer=TRUE)
    }
    
    for(tu in uniqueExon){
      ind <- which(rawexp[ ,"tu"] == tu)
      if(env == 1) axis(3, at=mean(ind), labels=tu, mgp=c(2.25, 0.5, 0), tick=FALSE, line=0)
      
      if(tu == uniqueExon[length(uniqueExon)]){
        rect(min(ind)-0.5, -3, fnPart[ ,2]+0.5, 20, density = 15, angle = 30, col = "azure2", border = "azure3")
        rect(fnPart[ ,2]+0.5, -3, max(ind)+0.5, 20, density = 15, angle = 150, col = "azure2", border = "azure3")
        
        #use median to test AS in 5'|3' site
        lines(c(fnPart[ ,2]+0.5, max(ind)+0.5), rep(median(unlist(newexp[exonID[exonID %in% (fnPart[ ,2]+1):max(ind)], ind_env])), 2), col=1, lty=2, lwd=1.5)
        
        if(fnPart[ ,env+2] >= asThre){
          text(x=(min(ind)+fnPart[ ,2])/2, y=max(newexp)+0.1, labels=round(fnPart[ ,env+2], digits=1), col=env+1, cex=1)
          lines(c(min(ind)-0.5, fnPart[ ,2]+0.5), rep(median(unlist(newexp[exonID[exonID %in% ind[1]:fnPart[ ,2]], ind_env])), 2), col=env+1, lwd=2)
        } else{
          if(round(fnPart[ ,env+2], digits=2) == asThre) text(x=(min(ind)+fnPart[ ,2])/2, y=max(newexp)+0.1, labels=round(fnPart[ ,env+2], digits=1), col=1, cex=1)
          lines(c(min(ind)-0.5, fnPart[ ,2]+0.5), rep(median(unlist(newexp[exonID[exonID %in% ind[1]:fnPart[ ,2]], ind_env])), 2), col=env+1, lty=2, lwd=1.5)
        }
      } else if(any(exonID %in% ind) && median(unlist(newexp[exonID[exonID %in% ind], ])) >= 5) lines(c(min(ind)-0.5, max(ind)+0.5), rep(median(unlist(newexp[exonID[exonID %in% ind], ind_env])), 2), col=1, lty=2, lwd=1.5)
    }
    
    for(p in 1:nrow(rawexp)){
      if(!p %in% ind_tu) rect((p-0.5), -3, (p+0.5), max(newexp[ ,ind_env])*2, col=grey(0.85), border="transparent")#background for introns
      if(p %in% probes_dir) points(rep(p, length(ind_env))+runif(length(ind_env), min=-0.05, max=0.05), newexp[p, ind_env], t='p', pch=20, cex=0.75)
    }
    box()
  }
  axis(1, at=probes_dir, labels=row.names(rawexp)[probes_dir], mgp=c(2.25, 0.5, 0), cex.axis=1, las=2, tck=0.02)
}


#plot most sig 3'AS
asMatrix <- NULL
for(chr in 1:5){
  if(file.exists(paste0("Data/AS/3'AS_chr", chr, "_wt_less_p2.txt"))){
    asMatrix <- rbind(asMatrix, read.table(paste0("Data/AS/3'AS_chr", chr, "_wt_less_p2.txt"), row.names=NULL))
  } else cat("chr", chr, "NO test!\n")
}

res <- NULL
for(e in 1:nrow(asMatrix)){
  res <- c(res, asMatrix[e, 3]*asMatrix[e, 4]*asMatrix[e, 5]*asMatrix[e, 6])
}
mostSig <- unique(asMatrix[order(res, decreasing=TRUE), 1])
for(filename in mostSig[1:5]){
  png(filename = paste0("Data/AS/plot/3'AS/", filename, "_3'AS_wt_mS.png"), width = 960, height = 1728, bg = "white")
  plot_3AS(filename, asThre=4.56)
  dev.off()
}
#plot around thre 3'AS
asMatrix <- NULL
for(chr in 1:5){
  if(file.exists(paste0("Data/AS/3'AS_chr", chr, "_wt_less_p2.txt"))){
    asMatrix <- rbind(asMatrix, read.table(paste0("Data/AS/3'AS_chr", chr, "_wt_less_p2.txt"), row.names=NULL))
  } else cat("chr", chr, "NO test!\n")
}

asThre=4.56
res <- NULL
for(e in 1:nrow(asMatrix)){
  res <- c(res, any(round(asMatrix[e, 3:6], digits=2) > (asThre-0.1) && round(asMatrix[e, 3:6], digits=2) < (asThre+0.1)))
}
aroundThre <- unique(asMatrix[res, 1])
for(filename in aroundThre[1:5]){
  png(filename = paste0("Data/AS/plot/3'AS/", filename, "_3'AS_wt_aT.png"), width = 960, height = 1728, bg = "white")
  plot_3AS(filename, asThre=4.56)
  dev.off()
}


#plot all 3'AS genes
for(chr in 1:5){
  load(file=paste0("Data/AS/3'AS_chr", chr, "_wt_p2.Rdata"))
  genenames <- seGeneList$mixEnv
  
  for(filename in genenames){
    png(filename = paste0("Data/AS/plot/3'AS/", filename, "_3'AS.png"), width = 960, height = 1728, bg = "white")
    plot_3AS(filename, asThre=4.56)
    dev.off()
  }
}



#*********************************************************** find best examples ***********************************************************#
getOne <- function(aa, cond = 1, cutoff=30){
  bb <- aa[which(aa[ ,cond] > cutoff), ]
  l10 <- sort(apply(bb[ ,(1:4)[-cond]], 1, sum), index.return=T)$ix[1:10]#sort(apply(bb[ ,-cond], 1, sum), index.return=T)$ix[1:10]
  bb[l10, ]
}

getOne <- function(aa, cond = 1, cutoff=30){
  bb <- aa[which(aa[ ,cond] > cutoff), ]
  l10 <- NULL
  for(r in 1:nrow(bb)){
    if(any(bb[r, -cond] <= (cutoff - 15))) l10 <- c(l10, r)
  }
  bb[l10, ]
}

for(chr in 1:5){
  #psmatrix <- read.table(paste0("Data/53terminalAS/53terminalAS_chr", chr, "_ttest.txt"), row.names=1, header=T)
  psmatrix <- read.table(paste0("Data/53terminalAS/53terminalAS_chr", chr, "_wtest.txt"), row.names=1, header=T)
  
  posPerfectEg <- NULL
  for(env in 1:4){
    posPerfectEg <- c(posPerfectEg, unlist(lapply(strsplit(rownames(getOne(psmatrix, env, 30)), "_"), "[[", 1)))
  }
  cat("I'm chr", chr, ", I have\n", posPerfectEg, "\n")
  for(filename in sort(unique(posPerfectEg))){
    #chr <- gsub("AT", "", unlist(lapply(strsplit(filename, "G"), "[[", 1)))
    plotcExonExp(chr, filename, ps_threshold = 5.03)
  }
}
