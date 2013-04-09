#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 8-04-2013
# first written: 8-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("D:/Arabidopsis Arrays")
#load environment file
env <- read.table("Data/ann_env.txt")
#load genomap
ann_m <- read.table("refined map/map.txt")


################################################################# panel 1 #################################################################
#list genes containing retention under any environment
retentionSum <- function(threshould=6, IRfile){
  retentionlist <- NULL
  for(env in 1:4){
    #retained <- the names of all retented introns of this chr
    retained <- rownames(IRfile[which(IRfile[ ,(env*4-1)] < threshould), ])
    #uniqueRetained <- the names of genes containing retented introns of this chr
    UniqueRetained <- unique(as.character(unlist(lapply(strsplit(retained,"_"),"[[",1))))
    
    retentionlist <- c(retentionlist, UniqueRetained)
  }
  retentionlist <- unique(retentionlist)
}

#choose genes randomly
geneID <- function(seed=1, ng=10, ...){
  set.seed(seed)
  
  retainedGenes <- retentionSum(threshould=6, IRfile)
  ind <- sample(1:length(retainedGenes), ng, replace=FALSE)
  orderind <- ind[order(ind, decreasing=FALSE)]
  filenames <- retainedGenes[orderind]
  
  return(filenames)
}

#select either good exons of right direction or introns of right direction
ProbesSelect <- function(chr, filename, exp_data=rawexp){
  #check the expression of exons of right direction
  load(paste("D:/Arabidopsis Arrays/Data/chr", chr, "_norm_hf_cor/", "Classification_chr", chr, "_norm_hf_cor.Rdata", sep=""))
  class <- res[[paste(filename, "_QTL.txt", sep="")]]#res[["AT1G01010_QTL.txt"]]
  
  probes <- NULL
  
  #decide which is right direction
  if(unique(exp_data[,"strand"]) == "sense"){
    direction_id <- which(exp_data[, "direction"] == "reverse")
  }
  if(unique(exp_data[,"strand"]) == "complement"){
    direction_id <- which(exp_data[, "direction"] == "forward")
  }
  
  #list good exons of right direction and introns of right direction
  probesDirExp <- direction_id[which(direction_id %in% class$goodP)]
  probes$In <- probesDirExp[which(probesDirExp %in% class$introP)]
  probes$Ex <- probesDirExp[-which(probesDirExp %in% class$introP)]
  
  return(probes)
}

#plot for expression intensity
makePlot_Exp <- function(filename, IRfile, rawexp, newexp, nprobes, env, ind_tu, ...){
  plot(c(0.5, nprobes+0.5), c(round(min(newexp)-0.5), round(max(newexp)+0.5)), xaxt='n', xlab="", ylab="Intensity", cex.axis=1, cex.lab=1.5, cex.main=2, las=1, mgp=c(3,1,0), tck=-0.017, t="n", main=filename)
  
  #retainedIn <- which(rawexp[,"tu"] %in% unique(as.character(unlist(lapply(strsplit(rownames(IRfile)[which(grepl(filename, rownames(IRfile)))],"_"),"[[",2)))))
  probesid <- unlist(ProbesSelect(chr, filename, exp_data=rawexp))
  
  for(p in 1:nprobes){
    #background for introns
    if(!p %in% ind_tu){
      rect((p-0.5), -3, (p+0.5), max(newexp)* 2.1, col=grey(0.85), border = "transparent")
    }
    #mark retained introns########################################## U N S U R E !!!!!*******************************************************
    #if(p %in% retainedIn){#####Here!!!!!
      #rect((p-0.5), -3, (p+0.5), max(newexp)* 2.1, col=grey(0.75), border = "transparent")
    #}
    #points of expression of probes used to check intron retention
    if(p %in% probesid){
      points(rep(p,148)+0.06*as.numeric(env[,2])-0.15, newexp[p,], t='p', col=env[,2], pch=20, cex=0.7)
    }
  }
  
  box();
}


################################################################# panel 2 #################################################################
#separate markers by chr
getProbesOnChr <- function(ann_m, chrs = 1){
  which(ann_m[,1] == chrs)
}

#plot for eQTL
makePlot_eQTL <- function(newexp, qtl, nprobes, ind_tu, lodThreshold = 5, lchr){
  plot(c(0.5, nprobes+0.5),c(0.5, 5.5), xaxt='n', xlab="", yaxt='n', ylab="eQTL", cex.axis=1, cex.lab=1.5, las=1, mgp=c(3,1,0), t="n")
  axis(2, at=1:5, labels=paste("chr", c("I", "II", "III", "IV", "V"),sep=""), cex.axis=1, las=2, tck=-0.035, mgp=c(2,0.5,0))
  
  #background for introns
  for(p in 1:nprobes){
    if(!p %in% ind_tu){
      rect((p-0.5), -3, (p+0.5), max(newexp)* 2.1, col=grey(0.85), border = "transparent")
    }
  }
  
  #line out which chr current gene are on
  abline(h=lchr, col="burlywood3", lty=8, lwd=3)
  points(-1, lchr, pch=">", cex=5, col="burlywood3")
  
  #rectangles for eQTL >= cutoff
  for(p in 1:nprobes){
    for(chrs in 1:5){
      if(any(qtl[p, getProbesOnChr(ann_m, chrs)] >= lodThreshold)){
        #cat("->",(round(max(qtl[p, getProbesOnChr(ann_m, chrs)]),d=0)+1) - lodThreshold,"\n")
        #points(p, chrs, pch=15,cex=1)
        rect((p-0.5),(chrs-0.4),(p+0.5),(chrs+0.4), col='tan2', border = "sienna3")
      }
    }
  }
  
  box();
}


################################################################# panel 3 #################################################################
#means of individuals under each environment for every probe
envTtest2 <- function(newexp, env, p){
  env1 <- mean(as.numeric(newexp[p, which(as.numeric(env[,2])==1)]))
  env2 <- mean(as.numeric(newexp[p, which(as.numeric(env[,2])==2)]))
  env3 <- mean(as.numeric(newexp[p, which(as.numeric(env[,2])==3)]))
  env4 <- mean(as.numeric(newexp[p, which(as.numeric(env[,2])==4)]))
  #lgp <- -log10(c(t.test(env1, mu=mean(c(env2,env3,env4)))$p.value,t.test(env2, mu=mean(c(env1,env3,env4)))$p.value,t.test(env3, mu=mean(c(env1,env2,env4)))$p.value,t.test(env4, mu=mean(c(env1,env2,env3)))$p.value))
  #env_mean <- mean(c(env1,env2,env3,env4))
  #lgp <- -log10(c(t.test(env1, mu=env_mean)$p.value,t.test(env2, mu=env_mean)$p.value,t.test(env3, mu=env_mean)$p.value,t.test(env4, mu=env_mean)$p.value))
  return(c(env1,env2,env3,env4))
}

#load color package and design color pattern
library("colorspace", lib.loc="C:/R/win-library/2.15")
cols <- diverge_hcl(31, h=c(195,330), c = 95, l = c(20, 90), power = 1.25)

#plot of environment comparison
makePlot_Env <- function(newexp, nprobes){
  plot(c(0.5, nprobes+0.5),c(0.5,4.5), xaxt='n', xlab="", yaxt='n', ylab="Env", cex.axis=1, cex.lab=1.5, las=1, mgp=c(3,1,0), t="n")
  axis(1, at=1:nprobes, labels=row.names(newexp), cex.axis=1, las=2, tck=0.035, mgp=c(2,1,0))
  axis(2, at=1, labels="6H", cex.axis=1, col.axis='black', las=2, tck=-0.035, mgp=c(2,0.5,0))
  axis(2, at=2, labels="DA", cex.axis=1, col.axis='red', las=2, tck=-0.035, mgp=c(2,0.5,0))
  axis(2, at=3, labels="DF", cex.axis=1, col.axis='green', las=2, tck=-0.035, mgp=c(2,0.5,0))
  axis(2, at=4, labels="RP", cex.axis=1, col.axis='blue', las=2, tck=-0.035, mgp=c(2,0.5,0))
  for(p in 1:nprobes){
    #pForCol <- round(envTtest(newexp,env,p)-0.5)+1
    pForCol <- ((envTtest2(newexp, env, p) - mean(as.numeric(newexp[p,])))*7.5) + 16
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
singlePlot <- function(rawexp, qtl, filename, ...){
  newexp <- rawexp[,17:164]
  ind_tu <- grep("tu", rawexp[,"tu"])
  nprobes <- nrow(qtl)
  st <- proc.time()
  par(fig=c(0,1,0.5,1), mar=c(0,5,5,2), oma=c(0,0,0,0.5))
  makePlot_Exp(filename, IRfile, rawexp, newexp, nprobes, env, ind_tu)
  par(fig=c(0,1,0.3,0.5), mar=c(0,5,0,2), oma=c(0,0,0,0.5), new=T)
  makePlot_eQTL(newexp, qtl, nprobes, ind_tu, lodThreshold = 5, lchr=chr)
  par(fig=c(0,1,0,0.3), mar=c(5,5,0,2), oma=c(0,0,0,0.5), new=T)
  makePlot_Env(newexp, nprobes)
  et <- proc.time()
  cat("Done with", filename, "after:",(et-st)[3],"secs\n")
}
#singlePlot(rawexp, qtl, plotname="Plot", lodThreshold = 4, lchr)


#makePlot(location = "C:/Arabidopsis Arrays/Data/chr1_norm_hf_cor/", lodThreshold = 4, lchr=1)
#lodThreshold = 4 <- -log10(0.05/716)=4.146128





makePlot <- function(chr, seed=1, ng=10, threshould=6){
  location = paste0("D:/Arabidopsis Arrays/Data/chr", chr, "_norm_hf_cor/")
  IRfile <- read.table(paste("Data/intronRetention/intronRetention_chr", chr, ".txt", sep=""), row.name=1, header=T)
  id <- geneID(seed=1, ng=10)
  for(filename in id){ #Please note: filename = AT1G01010!!!
    if(!file.exists(paste("Data/intronRetention/", filename, ".png", sep=""))){
      rawexp <- read.table(paste(location, filename, ".txt", sep=""), header=TRUE, row.names=1)
      qtl <- read.table(paste(location, filename, "_QTL.txt", sep=""), header=TRUE, row.names=1)
      
      png(file = paste("Data/intronRetention/", filename, ".png", sep=""), bg="white", width=1024, height=1024)
      singlePlot(rawexp, qtl, filename)
      dev.off()
    }
    else{
      cat("Skipping", filename," because it exists\n")
    }
  }
}

########################## chr unchanged!!!
for(chr in 1:5){
  makePlot(chr, seed=1, ng=10, threshould=6)
}
