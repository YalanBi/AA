#P3-Generate pictures 3in1

setwd("C:/Arabidopsis Arrays")
geno <- read.table("Data/genotypes_n.txt", row.names=1)
ann <- read.table("annotation_sample.txt", colClasses="character")
menvironment <- ann[9:172,3]
map.fast <- function(x, geno, pheno, menvironment){
	res <- NULL
	models  <- aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno[,x]))
	modelinfo <- summary(models)
	res$env <- unlist(lapply(modelinfo,"[",1,5),use.names=T)
	res$qtl <- unlist(lapply(modelinfo,"[",2,5),use.names=T)
	res
}



rr=c(244,239,225,215,209,193,181,166,151,130)
gg=c(228,204,174,160,146,117,94,58,44,45)
bb=c(176,140,109,105,102,91,84,74,70,68)

mycolor=NULL
for ( i in 1:10){
	mycolor=c(mycolor,rgb(rr[i],gg[i],bb[i],maxColorValue=255)  )
}
mycolor



st  <- proc.time()
for(filename in dir("Data/gene_data")[grepl(".txt",dir("Data/gene_data"))]){
  cat("Loading", filename,"\n")
  Pheno <- read.table(paste("Data/gene_data/",filename,sep=""))
  pheno <- t(Pheno[,25:188])
  
  ind_tu <- grep("tu", Pheno[,9])
  
  resmatrix_g <- NULL
  for(x in 1:ncol(geno)){
    resmatrix_g <- cbind(resmatrix_g, -log(map.fast(x, geno, pheno, menvironment)$qtl))
  }
  #write.table(resmatrix, file=paste("Data/gene_eff~lm", gsub(".txt","_G.txt",filename),sep=""), row.names = colnames(pheno), col.names = paste("g", 1:69, sep=""))

  #generate a matrix of the means of environments
  res <- NULL
  mean_env <- NULL
  for(x in 1:ncol(pheno)){
    mean_env <- c(mean(pheno[which(as.character(menvironment) == "6H"),x]),mean(pheno[which(as.character(menvironment) == "Dry_Fresh"),x]),mean(pheno[which(as.character(menvironment) == "Dry_AR"),x]),mean(pheno[which(as.character(menvironment) == "RP"),x]))
    mean_env <- c(mean_env, mean(mean_env), sd(mean_env), (mean_env-mean(mean_env))/sd(mean_env))
    res <- cbind(res,mean_env)
  }
  colnames(res) <- colnames(pheno)
  rownames(res) <- c("6H", "Dry_Fresh", "Dry_AR", "RP", "gmean", "gsd", "(1-m)/sd",  "(2-m)/sd", "(3-m)/sd", "(4-m)/sd")
  
  if(ncol(pheno) >= 4){
    #png(file = paste("Data/Pic/",gsub(".txt","_p.png",filename),sep=""), bg="white", width=1024, height=1024)
    par(fig=c(0,1,0.7,1), mar=c(0,2.5,2,0), oma=c(0,0,0,0.5))
    plot(c(1,ncol(pheno)),c(min(pheno) * .9,max(pheno)* 1.1), xaxt='n', xlab="", ylab="Intensity", cex.axis=0.75, cex.lab=0.9, cex.main=1, las=1, mgp=c(1.5,0.5,0), t="n", main=gsub("_G.txt","",filename))
    for(p in 1:ncol(pheno)){
      if(p %in% ind_tu){
        points(rep(p,164)+((runif(164)-0.5)/2), pheno[,p], t='p', col='blue',pch=20,cex=0.7)
      } else {
        rect((p-0.5), 0,(p+0.5),max(pheno)* 2.1,col=grey(0.85),border = "transparent")
		points(rep(p,164)+((runif(164)-0.5)/2.5), pheno[,p], t='p', col='blue',pch=20,cex=0.7)
      }
    }
	s <- 1
	for(x in 1:ncol(pheno)){
	if(Pheno[x,"tu"] != Pheno[s,"tu"]){
		lines(c(s-0.5,x-0.5),c(mean(pheno[,s:x-1]),mean(pheno[,s:x-1])))
		s <- x
		}
	}
	lines(c(s-0.5,x+0.5),c(mean(pheno[,s:x]),mean(pheno[,s:x])))
	box();
    
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
	
    par(fig=c(0,1,0,0.2), mar=c(2.5,2.5,0,0), oma=c(0,0,0,0.5), new=T)
    plot(c(1,ncol(pheno)),c(0.5,4.5), xlab="Probes", ylab="Env", cex.axis=0.75, cex.lab=0.9, las=1, mgp=c(1.5,0.5,0), t="n")
    for(a in 1:ncol(pheno)){
      gmean <- mean(res[,a])
      gsd   <- sd(res[,a])
      rect((a-0.5),0.6,(a+0.5),1.4,col=mycolor[round(res[7,a]*3+5.5)],border = "transparent")
      rect((a-0.5),1.6,(a+0.5),2.4,col=mycolor[round(res[8,a]*3+5.5)],border = "transparent")
      rect((a-0.5),2.6,(a+0.5),3.4,col=mycolor[round(res[9,a]*3+5.5)],border = "transparent")
      rect((a-0.5),3.6,(a+0.5),4.4,col=mycolor[round(res[10,a]*3+5.5)],border = "transparent")
    }
    #dev.off()
  }
  cat("Done in:",(proc.time()-st)[3],"seconds\n")
  #return()
}

points(1,5,pch=19,col=rgb(.25,.75,.75))
points(2,5,pch=19,col=rgb(1,.5,1))
points(3,5,pch=19,col=rgb(0.5,.75,.25))
points(4,5,pch=19,col=rgb(1,.5,.25))








#means for all genes
envall<- t(read.table("exp_ann30K.txt")[,25:188])
resall <- NULL
mean_envall <- NULL
for(p_a in 1:ncol(envall)){
  mean_envall <- c(mean(envall[which(as.character(menvironment) == "6H"),p_a]),mean(envall[which(as.character(menvironment) == "Dry_Fresh"),p_a]),mean(envall[which(as.character(menvironment) == "Dry_AR"),p_a]),mean(envall[which(as.character(menvironment) == "RP"),p_a]))
  mean_envall <- c(mean_envall, mean(mean_envall), sd(mean_envall), (mean_envall-mean(mean_envall))/sd(mean_envall))
  resall <- rbind(resall,mean_envall)
  #return()
}
rownames(resall) <- colnames(envall)
colnames(resall) <- c("6H", "Dry_Fresh", "Dry_AR", "RP", "gmean", "gsd", "(1-m)/sd",  "(2-m)/sd", "(3-m)/sd", "(4-m)/sd")
write.table(resall,file="means_evn_all_2.txt")