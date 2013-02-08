#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 29-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("C:/Arabidopsis Arrays")
env <- read.table("Data/ann_env.txt")
ann_m <- read.table("refined map/map.txt")

makePlot_Exp <- function(fn_exp, rawexp, newexp, nprobes, env, ind_tu){
  plot(c(0.5, nprobes+0.5), c(round(min(newexp)-0.5), round(max(newexp)+0.5)), xaxt='n', xlab="", ylab="Intensity", cex.axis=1, cex.lab=1.5, cex.main=2, las=1, mgp=c(3,1,0), tck=-0.017, t="n", main=gsub(".txt", "", fn_exp))
  s <- 1
  for(p in 1:nprobes){
    if(!p %in% ind_tu){
      rect((p-0.5), -3, (p+0.5), max(newexp)* 2.1, col=grey(0.85), border = "transparent")
    }
    points(rep(p,148)+0.06*as.numeric(env[,2])-0.15, newexp[p,], t='p', col=env[,2], pch=20, cex=0.7)
    if(rawexp[p,"tu"] != rawexp[s,"tu"]){
      lines(c(s-0.5, p-0.5), c(mean(as.matrix(newexp[s:(p-1),])), mean(as.matrix(newexp[s:(p-1),]))))
      s <- p
    }
  }
  lines(c(s-0.5, p+0.5), c(mean(as.matrix(newexp[s:p,])), mean(as.matrix(newexp[s:p,]))))
  box();
}

getProbesOnChr <- function(ann_m, chr = 1){
  which(ann_m[,1] == chr)
}

makePlot_eQTL <- function(newexp, qtl, nprobes, ind_tu, lodThreshold = 5, chrs = 1:5, lchr){
  plot(c(0.5, nprobes+0.5),c(0.5, 5.5), xaxt='n', xlab="", yaxt='n', ylab="eQTL", cex.axis=1, cex.lab=1.5, las=1, mgp=c(3,1,0), t="n")
  axis(2, at=1:5, labels=paste("chr", c("I", "II", "III", "IV", "V"),sep=""), cex.axis=1, las=2, tck=-0.035, mgp=c(2,0.5,0))
  for(p in 1:nprobes){
    if(!p %in% ind_tu){
      rect((p-0.5), -3, (p+0.5), max(newexp)* 2.1, col=grey(0.85), border = "transparent")
    }
  }
  abline(h=lchr, col="burlywood3", lty=8, lwd=3)
  points(-1,lchr, pch=">", cex=5, col="burlywood3")
  for(p in 1:nprobes){
    for(chr in chrs){
      if(any(qtl[p, getProbesOnChr(ann_m, chr)] >= lodThreshold)){
        #cat("->",(round(max(qtl[p, getProbesOnChr(ann_m, chr)]),d=0)+1) - lodThreshold,"\n")
        #points(p, chr, pch=15,cex=1)
        rect((p-0.5),(chr-0.4),(p+0.5),(chr+0.4), col='tan2', border = "sienna3")
      }
    }
  }
  
  box();
}

envTtest <- function(newexp,env, p){
  env1 <- as.matrix(newexp[p,which(as.numeric(env[,2])==1)])
  env2 <- as.matrix(newexp[p,which(as.numeric(env[,2])==2)])
  env3 <- as.matrix(newexp[p,which(as.numeric(env[,2])==3)])
  env4 <- as.matrix(newexp[p,which(as.numeric(env[,2])==4)])
  lgp <- -log10(c(t.test(env1, mu=mean(c(env2,env3,env4)))$p.value,t.test(env2, mu=mean(c(env1,env3,env4)))$p.value,t.test(env3, mu=mean(c(env1,env2,env4)))$p.value,t.test(env4, mu=mean(c(env1,env2,env3)))$p.value))
  #env_mean <- mean(c(env1,env2,env3,env4))
  #lgp <- -log10(c(t.test(env1, mu=env_mean)$p.value,t.test(env2, mu=env_mean)$p.value,t.test(env3, mu=env_mean)$p.value,t.test(env4, mu=env_mean)$p.value))
  return(lgp)
}

envTtest2 <- function(newexp,env, p){
  env1 <- mean(as.numeric(newexp[p,which(as.numeric(env[,2])==1)]))
  env2 <- mean(as.numeric(newexp[p,which(as.numeric(env[,2])==2)]))
  env3 <- mean(as.numeric(newexp[p,which(as.numeric(env[,2])==3)]))
  env4 <- mean(as.numeric(newexp[p,which(as.numeric(env[,2])==4)]))
 # lgp <- -log10(c(t.test(env1, mu=mean(c(env2,env3,env4)))$p.value,t.test(env2, mu=mean(c(env1,env3,env4)))$p.value,t.test(env3, mu=mean(c(env1,env2,env4)))$p.value,t.test(env4, mu=mean(c(env1,env2,env3)))$p.value))
  #env_mean <- mean(c(env1,env2,env3,env4))
  #lgp <- -log10(c(t.test(env1, mu=env_mean)$p.value,t.test(env2, mu=env_mean)$p.value,t.test(env3, mu=env_mean)$p.value,t.test(env4, mu=env_mean)$p.value))
  return(c(env1,env2,env3,env4))
}

mycolor <- function(){
  rr=c(244,239,225,215,209,193,181,166,151,130)
  gg=c(228,204,174,160,146,117,94,58,44,45)
  bb=c(176,140,109,105,102,91,84,74,70,68)
  collist <- NULL
  for (i in 1:10){
    collist <- c(collist,rgb(rr[i],gg[i],bb[i],maxColorValue=255))
  }
  return(collist)
}

mycolz <- function(){
  colz <- NULL
  for(x in 1:15){ colz <- c(colz, rgb(0.9,x/20,0.9)) }
  colz <- c(colz, rgb(0.9,0.9,0.9))
  for(x in 1:15){ colz <- c(colz, rgb((20-x)/20,0.9,0.9)) }
  return(colz)
}

cols <- diverge_hcl(31, h=c(195,330), c = 95, l = c(20, 90), power = 1.25)

makePlot_Env <- function(newexp, env, nprobes){
  plot(c(0.5, nprobes+0.5),c(0.5,4.5), xaxt='n', xlab="", yaxt='n', ylab="Env", cex.axis=1, cex.lab=1.5, las=1, mgp=c(3,1,0), t="n")
  axis(1, at=1:nprobes, labels=row.names(newexp), cex.axis=1, las=2, tck=0.035, mgp=c(2,1,0))
  axis(2, at=1, labels="6H", cex.axis=1, col.axis='black', las=2, tck=-0.035, mgp=c(2,0.5,0))
  axis(2, at=2, labels="DA", cex.axis=1, col.axis='red', las=2, tck=-0.035, mgp=c(2,0.5,0))
  axis(2, at=3, labels="DF", cex.axis=1, col.axis='green', las=2, tck=-0.035, mgp=c(2,0.5,0))
  axis(2, at=4, labels="RP", cex.axis=1, col.axis='blue', las=2, tck=-0.035, mgp=c(2,0.5,0))
  #oamean <- mean(unlist(newexp))
  for(p in 1:nprobes){
    #pForCol <- round(envTtest(newexp,env,p)-0.5)+1
    pForCol <- ((envTtest2(newexp,env,p) - mean(as.numeric(newexp[p,])))*7.5) + 16
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

#Please note: Filename = QTL file !!
singlePlot <- function(rawexp, qtl, fn_exp="Plot", lodThreshold = 5, lchr){
  newexp <- rawexp[,17:164]
  ind_tu <- grep("tu", rawexp[,"tu"])
  nprobes <- nrow(qtl)
  st <- proc.time()
  par(fig=c(0,1,0.5,1), mar=c(0,5,5,2), oma=c(0,0,0,0.5))
  makePlot_Exp(fn_exp, rawexp, newexp, nprobes, env, ind_tu)
  par(fig=c(0,1,0.3,0.5), mar=c(0,5,0,2), oma=c(0,0,0,0.5), new=T)
  makePlot_eQTL(newexp, qtl, nprobes, ind_tu, lodThreshold = lodThreshold, chrs = 1:5, lchr)
  par(fig=c(0,1,0,0.3), mar=c(5,5,0,2), oma=c(0,0,0,0.5), new=T)
  makePlot_Env(newexp, env, nprobes)
  et <- proc.time()
  cat("Done with plot after:",(et-st)[3],"secs\n")
}
#singlePlot(rawexp, qtl, fn_exp="Plot", lodThreshold = 4, lchr)

makePlot <- function(location = "C:/Arabidopsis Arrays/Data/chr1_norm_hf_cor/", ...){
  for(filename in dir(location)[grepl("_QTL",dir(location))]){ #Please note: Filename = QTL file !!
    fn_qtl <- filename
    fn_exp <- gsub("_QTL.txt",".txt", filename)
    fn_png <- gsub("_QTL.txt",".png", filename)
    if(!file.exists(paste(location, fn_png, sep=""))){
      rawexp <- read.table(paste(location, fn_exp, sep=""), header=TRUE, row.names=1)
      qtl <- read.table(paste(location, fn_qtl, sep=""), header=TRUE, row.names=1)
      png(file = paste(location, fn_png, sep=""), bg="white", width=1024, height=1024)
      singlePlot(rawexp, qtl, fn_exp, ...)
      dev.off()
    }else{
      cat("Skipping", fn_qtl," because it exists\n")
    }
  }
}

makePlot(location = "C:/Arabidopsis Arrays/Data/chr1_norm_hf_cor/", lodThreshold = 4, lchr=1)
#lodThreshold = 4 <- -log10(0.05/716)=4.146128
