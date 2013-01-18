#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 16-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


geno <- read.table("Data/genotypes_n.txt", row.names=1)
ann <- read.table("annotation_sample.txt", colClasses="character")
menvironment <- ann[9:172,3]


makePlot <- function(filename, x, ...){
	plot(x ,...)
}


    #image(t(resmatrix_g))
    par(fig=c(0,1,0.2,0.7), mar=c(0,2.5,0,0), oma=c(0,0,0,0.5), new=T)
    plot(c(1,ncol(pheno)),c(1,ncol(resmatrix_g)), xaxt='n', xlab="", ylab="eQTL", cex.axis=0.75, cex.lab=0.9, las=1, mgp=c(1.5,0.5,0), t="n")
    for(p in 1:nrow(resmatrix_g)){
		if(p %in% ind_tu){}
		else {
        rect((p-0.5),-2,(p+0.5),ncol(resmatrix_g)+5,col=grey(0.85),border = "transparent")
      }
      for(m in 1:ncol(resmatrix_g)){
        if(resmatrix_g[p,m] >= 4){
          points(p,m, pch=20,cex=1, col=rgb(0.4,0,0.6,0.5:(max(resmatrix_g)+1)/(max(resmatrix_g)+1))[(round(resmatrix_g[p,m])+1)])
        }
      }
	  
    }
    box();

points(1,5,pch=19,col=rgb(.25,.75,.75))
points(2,5,pch=19,col=rgb(1,.5,1))
points(3,5,pch=19,col=rgb(0.5,.75,.25))
points(4,5,pch=19,col=rgb(1,.5,.25))




setwd("C:/Arabidopsis Arrays")
rawexp <- read.table("Data/chr1/AT1G01010.txt", header=TRUE, row.names=1)
newexp <- rawexp[,17:164]
env <- read.table("Data/ann_env.txt")
ind_tu <- grep("tu", rawexp[,"tu"])


makePlot_Exp <- function(filename, rawexp, newexp, env, ind_tu, ...){
  s <- 1
  for(p in 1:nrow(newexp)){
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
par(fig=c(0,1,0.5,1), mar=c(0,2.5,2,0), oma=c(0,0,0,0.5))
plot(c(1, nrow(newexp)), c(min(newexp) * .9, max(newexp)* 1.1), xaxt='n', xlab="", ylab="Intensity", cex.axis=0.75, cex.lab=0.9, cex.main=1, las=1, mgp=c(1.5,0.5,0), t="n", main=filename)
makePlot_Exp(filename, rawexp, newexp, env, ind_tu)



ann_m <- read.table("refined map/map.txt")
getProbesOnChr <- function(ann_m, chr = 1){
  which(ann_m[,1] == chr)
}

makePlot_eQTL <- function(filename, newexp, qtl, ind_tu, lodThreshold = 4, chrs = 1:5){
  qtl <- read.table(paste("Data/chr1/", gsub(".txt", "_QTL.txt", filename), sep=""))
  nprobes <- nrow(qtl)
  for(chr in chrs){
    for(p in 1:nprobes){
      if(!p %in% ind_tu){
        rect((p-0.5), -3, (p+0.5), max(newexp)* 2.1, col=grey(0.85), border = "transparent")
      }
    }
  }
  for(chr in chrs){
    for(p in 1:nprobes){
      if(any(qtl[p, getProbesOnChr(ann_m, chr)] >= lodThreshold)){
        points(p,chr, pch=20,cex=1)
      }
    }
  }  
  box();
}
par(fig=c(0,1,0.2,0.5), mar=c(0,2.5,0,0), oma=c(0,0,0,0.5), new=T)
plot(c(1, nprobes),c(1,5), xaxt='n', xlab="", ylab="eQTL", cex.axis=0.75, cex.lab=0.9, las=1, mgp=c(1.5,0.5,0), t="n")
makePlot_eQTL(filename, newexp, qtl, ind_tu, lodThreshold = 4, chrs = 1:5)








envTtest <- function(newexp,env, p){
  env1 <- as.matrix(newexp[p,which(as.numeric(env[,2])==1)])
  env2 <- as.matrix(newexp[p,which(as.numeric(env[,2])==2)])
  env3 <- as.matrix(newexp[p,which(as.numeric(env[,2])==3)])
  env4 <- as.matrix(newexp[p,which(as.numeric(env[,2])==4)])
  env_mean <- mean(c(env1,env2,env3,env4))
  lgp <- -log(c(t.test(env1, mu=env_mean)$p.value,t.test(env2, mu=env_mean)$p.value,t.test(env3, mu=env_mean)$p.value,t.test(env4, mu=env_mean)$p.value))
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

makePlot_Env <- function(filename, newexp, env, p, ...){
  for(p in 1:nrow(newexp)){
    pForCol <- round(envTtest(newexp,env,p))+1
    rect((p-0.5),0.6,(p+0.5),1.4,col=mycolor()[pForCol[1]],border = "transparent")
    rect((p-0.5),1.6,(p+0.5),2.4,col=mycolor()[pForCol[2]],border = "transparent")
    rect((p-0.5),2.6,(p+0.5),3.4,col=mycolor()[pForCol[3]],border = "transparent")
    rect((p-0.5),3.6,(p+0.5),4.4,col=mycolor()[pForCol[4]],border = "transparent")
  }
}
par(fig=c(0,1,0,0.2), mar=c(2.5,2.5,0,0), oma=c(0,0,0,0.5), new=T)
plot(c(1,nrow(newexp)),c(0.5,4.5), xlab="Probes", ylab="Env", cex.axis=0.75, cex.lab=0.9, las=1, mgp=c(1.5,0.5,0), t="n")
makePlot_Env(filename, newexp, env, p)
