#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 24-01-2013
# first written: 22-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("C:/Arabidopsis Arrays")
env <- read.table("Data/ann_env.txt")
ann_m <- read.table("refined map/map.txt")

makePlot_Exp <- function(fn_exp, rawexp, newexp, nprobes, env, ind_tu){
  plot(c(0.5, nprobes+0.5), c(min(newexp) * .9, max(newexp)* 1.1), xaxt='n', xlab="", ylab="Intensity", cex.axis=0.75, cex.lab=1.2, cex.main=1.5, las=1, mgp=c(1.5,0.5,0), t="n", main=fn_exp)
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

makePlot_eQTL <- function(newexp, qtl, nprobes, ind_tu, lodThreshold = 5, chrs = 1:5){
  plot(c(0.5, nprobes+0.5),c(0.5, 5.5), xaxt='n', xlab="", ylab="eQTL", cex.axis=0.75, cex.lab=1.2, las=1, mgp=c(1.5,0.5,0), t="n")
  for(p in 1:nprobes){
    if(!p %in% ind_tu){
      rect((p-0.5), -3, (p+0.5), max(newexp)* 2.1, col=grey(0.85), border = "transparent")
    }
    for(chr in chrs){
      if(any(qtl[p, getProbesOnChr(ann_m, chr)] >= lodThreshold)){
        #cat("->",(round(max(qtl[p, getProbesOnChr(ann_m, chr)]),d=0)+1) - lodThreshold,"\n")
        points(p, chr, pch=15,cex=1)
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
  env_mean <- mean(c(env1,env2,env3,env4))
  lgp <- -log10(c(t.test(env1, mu=env_mean)$p.value,t.test(env2, mu=env_mean)$p.value,t.test(env3, mu=env_mean)$p.value,t.test(env4, mu=env_mean)$p.value))
  return(lgp)
}
mycolor <- function(){
  rr=c(244,239,225,215,209,193,181,166,151,130)
  gg=c(228,204,174,160,146,117,94,58,44,45)
  bb=c(176,140,109,105,102,91,84,74,70,68)
  collist <- NULL
  for ( i in 1:10){
    collist <- c(collist,rgb(rr[i],gg[i],bb[i],maxColorValue=255))
  }
  return(collist)
}

makePlot_Env <- function(newexp, env, nprobes){
  plot(c(0.5, nprobes+0.5),c(0.5,4.5), xlab="Probes", ylab="Env", cex.axis=0.75, cex.lab=1.2, las=1, mgp=c(1.5,0.5,0), t="n")
  for(p in 1:nprobes){
    pForCol <- round(envTtest(newexp,env,p)-0.5)+1
    pForCol[pForCol > 10] <- 10
    rect((p-0.5),0.6,(p+0.5),1.4,col=mycolor()[pForCol[1]],border = "transparent")
    rect((p-0.5),1.6,(p+0.5),2.4,col=mycolor()[pForCol[2]],border = "transparent")
    rect((p-0.5),2.6,(p+0.5),3.4,col=mycolor()[pForCol[3]],border = "transparent")
    rect((p-0.5),3.6,(p+0.5),4.4,col=mycolor()[pForCol[4]],border = "transparent")
  }
  box();
}

#Please note: Filename = QTL file !!
singlePlot <- function(rawexp, qtl, fn_exp="Plot", lodThreshold = 5){
  newexp <- rawexp[,16:163]
  ind_tu <- grep("tu", rawexp[,"tu"])
  nprobes <- nrow(qtl)
  st <- proc.time()
  par(fig=c(0,1,0.5,1), mar=c(0,2.5,2,0), oma=c(0,0,0,0.5))
  makePlot_Exp(fn_exp, rawexp, newexp, nprobes, env, ind_tu)
  par(fig=c(0,1,0.2,0.5), mar=c(0,2.5,0,0), oma=c(0,0,0,0.5), new=T)
  makePlot_eQTL(newexp, qtl, nprobes, ind_tu, lodThreshold = lodThreshold, chrs = 1:5)
  par(fig=c(0,1,0,0.2), mar=c(2.5,2.5,0,0), oma=c(0,0,0,0.5), new=T)
  makePlot_Env(newexp, env, nprobes)
  et <- proc.time()
  cat("Done with plot after:",(et-st)[3],"secs\n")
}
#singlePlot(rawexp, qtl, fn_exp="Plot", lodThreshold = 5)

makePlot <- function(location = "C:/Arabidopsis Arrays/Data/chr1/", ...){
  for(filename in dir(location)[grepl("_QTL",dir(location))]){ #Please note: Filename = QTL file !!
    fn_qtl <- filename
    fn_exp <- gsub("_QTL.txt",".txt", filename)
    fn_png <- gsub("_QTL.txt",".png", filename)
    if(!file.exists(paste(location, fn_png, sep=""))){
      rawexp <- read.table(paste(location, fn_exp, sep=""), header=TRUE, row.names=1)[res[[filename]]$gooddP,]
      qtl <- read.table(paste(location, fn_qtl, sep=""), header=TRUE, row.names=1)[res[[filename]]$gooddP,]
      png(file = paste(location, fn_png, sep=""), bg="white", width=1024, height=1024)
      singlePlot(rawexp, qtl, fn_exp, ...)
      dev.off()
    }else{
      cat("Skipping", fn_qtl," because it exists\n")
    }
  }
}

makePlot(location = "C:/Arabidopsis Arrays/Data/chr1/", lodThreshold = 5)
#lodThreshold = 5 <- -log10(0.01/716)
