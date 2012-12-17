#P3-Generate pictures 3in1

setwd("X:/Data/Desktop/Arabidopsis Arrays")
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
		res <- cbind(res,mean_env)
	}
	colnames(res) <- colnames(pheno)
	rownames(res) <- c("6H", "Dry_Fresh", "Dry_AR", "RP")

	if(ncol(pheno) >= 4){
		png(file = paste("Data/Pic/",gsub(".txt","_p.png",filename),sep=""), bg="white", width=1024, height=1024)
		op <- par(mfrow=c(3,1), mar=c(1,0,0,0), oma=c(2,2,0,0))
		plot(c(0,(ncol(pheno)+1)),c(min(pheno) * .9,max(pheno)* 1.1), xlab="", ylab="", t="n")
		for(p in 1:ncol(pheno)){
			if(p %in% ind_tu){
				points(rep(p,164)+((runif(164)-0.5)/2), pheno[,p], t='p', col='blue',pch=20,cex=0.7)
			} else {
				points(rep(p,164)+((runif(164)-0.5)/2), pheno[,p], t='p', col='gray',pch=20,cex=0.7)
			}
		}
		#image(t(resmatrix_g))
		plot(c(1,(nrow(resmatrix_g)+1)),c(1,(ncol(resmatrix_g)+1)), xlab="", ylab="", xlim = c(1,(nrow(resmatrix_g)+1)), ylim = c(0,(ncol(resmatrix_g)+1)), t="n")
		for(p in 1:nrow(resmatrix_g)){
			for(m in 1:ncol(resmatrix_g)){
				if(resmatrix_g[p,m] >= 4){
					points(p,m, pch=20, col=rgb(0.7,0,0:(max(resmatrix_g)+1)/(max(resmatrix_g)+1))[(round(resmatrix_g[p,m])+1)])
				}
			}
		}
		
		plot(c(0,(ncol(pheno)+1)),c(0.5,4.5), xlab="", ylab="", t="n")
		for(a in 1:ncol(pheno)){
			rect((a-0.5),0.5,(a+0.5),1.5,col=rgb(0:(max(res)+1)/(max(res)+1),1-0:(max(res)+1)/(max(res)+1),0)[(round(res[1,a])+1)],border = "transparent")
			rect((a-0.5),1.5,(a+0.5),2.5,col=rgb(0:(max(res)+1)/(max(res)+1),1-0:(max(res)+1)/(max(res)+1),0)[(round(res[2,a])+1)],border = "transparent")
			rect((a-0.5),2.5,(a+0.5),3.5,col=rgb(0:(max(res)+1)/(max(res)+1),1-0:(max(res)+1)/(max(res)+1),0)[(round(res[3,a])+1)],border = "transparent")
			rect((a-0.5),3.5,(a+0.5),4.5,col=rgb(0:(max(res)+1)/(max(res)+1),1-0:(max(res)+1)/(max(res)+1),0)[(round(res[4,a])+1)],border = "transparent")
		}
		dev.off()
	}
	cat("Done in:",(proc.time()-st)[3],"seconds\n")
	return()
}





makePlot <- function(filename, x, ...){
	plot(x ,...)
}
png(file = paste("Data/Pic/",gsub(".txt","_p.png",filename),sep=""), bg="white", width=1024, height=1024)
par(fig=c(0,1,0.7,1), mar=c(0,2.5,2,0), oma=c(0,0,0,0.5))
plot(c(1,ncol(pheno)),c(min(pheno) * .9,max(pheno)* 1.1), xaxt='n', xlab="", ylab="Intensity", cex.axis=0.5, cex.lab=0.75, cex.main=0.85, las=1, mgp=c(1.5,0.5,0), t="n", main=gsub("_G.txt","",filename))
for(p in 1:ncol(pheno)){
	if(p %in% ind_tu){
		points(rep(p,164)+((runif(164)-0.5)/2), pheno[,p], t='p', col='blue',pch=20,cex=0.7)
	} else {
		points(rep(p,164)+((runif(164)-0.5)/2), pheno[,p], t='p', col='gray',pch=20,cex=0.7)
	}
}
par(fig=c(0,1,0.2,0.7), mar=c(0,2.5,0,0), oma=c(0,0,0,0.5), new=T)
plot(c(1,ncol(pheno)),c(1,(ncol(resmatrix_g)+1)), xaxt='n', xlab="", ylab="eQTL", cex.axis=0.5, cex.lab=0.75, las=1, mgp=c(1.5,0.5,0), t="n")
for(p in 1:nrow(resmatrix_g)){
	for(m in 1:ncol(resmatrix_g)){
		if(resmatrix_g[p,m] >= 4){
			points(p,m, pch=20, col=rgb(0.5,0,0:(max(resmatrix_g)+1)/(max(resmatrix_g)+1))[(round(resmatrix_g[p,m])+1)])
		}#colors!
	}
}
par(fig=c(0,1,0,0.2), mar=c(2.5,2.5,0,0), oma=c(0,0,0,0.5), new=T)
plot(c(1,ncol(pheno)),c(0.5,4.5), xlab="Probes", ylab="Env", cex.axis=0.5, cex.lab=0.75, las=1, mgp=c(1.5,0.5,0), t="n")
for(a in 1:ncol(pheno)){
	rect((a-0.5),0.5,(a+0.5),1.5,col=rgb(0:(max(res)+1)/(max(res)+1),1-0:(max(res)+1)/(max(res)+1),0)[(round(res[1,a])+1)],border = "transparent")
	rect((a-0.5),1.5,(a+0.5),2.5,col=rgb(0:(max(res)+1)/(max(res)+1),1-0:(max(res)+1)/(max(res)+1),0)[(round(res[2,a])+1)],border = "transparent")
	rect((a-0.5),2.5,(a+0.5),3.5,col=rgb(0:(max(res)+1)/(max(res)+1),1-0:(max(res)+1)/(max(res)+1),0)[(round(res[3,a])+1)],border = "transparent")
	rect((a-0.5),3.5,(a+0.5),4.5,col=rgb(0:(max(res)+1)/(max(res)+1),1-0:(max(res)+1)/(max(res)+1),0)[(round(res[4,a])+1)],border = "transparent")
}





##############~~~~~~Here!!!~~~~~~#############
op <- par(mfrow=c(3,1), mar=c(0,0,0,0), oma=c(2,0,0,0))
plot(c(0,(ncol(pheno)+1)),c(1,10), xlab="", ylab="", t="n")
for(p in 1:nrow(pheno)){
	points(pheno[p,], t='p', col=Pheno[,9])
}
#image(t(resmatrix_g))
plot(c(1,(nrow(resmatrix_g)+1)),c(1,(ncol(resmatrix_g)+1)), xlab="", ylab="", t="n")
for(p in 1:nrow(resmatrix_g)){
	for(m in 1:ncol(resmatrix_g)){
		if(resmatrix_g[p,m] >= 4){
			points(p,m, pch=20, col=gray(0:max(resmatrix_g)/max(resmatrix_g))[(round(resmatrix_g[p,m])+1)])
		}
	}
}
#points(c(1:nrow(resmatrix_g)),rep(x,nrow(resmatrix_g)), pch=20, col=rgb(0:max(resmatrix_g)/max(resmatrix_g),1-0:max(resmatrix_g)/max(resmatrix_g),1)[(round(resmatrix_g[,x])+1)])
plot(c(0,(ncol(pheno)+1)),c(0.5,4.5), t="n")
for(a in 1:ncol(pheno)){
	rect((a-0.5),0.5,(a+0.5),1.5,col=rgb((0:(max(res)+1))/(max(res)+1),1-(0:(max(res)+1))/(max(res)+1),0.3)[(round(res[1,a])+1)],border = "transparent")
	rect((a-0.5),1.5,(a+0.5),2.5,col=rgb((0:(max(res)+1))/(max(res)+1),1-(0:(max(res)+1))/(max(res)+1),0.3)[(round(res[2,a])+1)],border = "transparent")
	rect((a-0.5),2.5,(a+0.5),3.5,col=rgb((0:(max(res)+1))/(max(res)+1),1-(0:(max(res)+1))/(max(res)+1),0.3)[(round(res[3,a])+1)],border = "transparent")
	rect((a-0.5),3.5,(a+0.5),4.5,col=rgb((0:(max(res)+1))/(max(res)+1),1-(0:(max(res)+1))/(max(res)+1),0.3)[(round(res[4,a])+1)],border = "transparent")
}





gray((0:(max(res)-min(res)+2))/(max(res)-min(res)+2))[(round(res[,a]-min(res))+1)]
col=rgb((0:(max(res)+1))/(max(res)+1),1-(0:(max(res)+1))/(max(res)+1),0.5)[(round(res[,a])+1)]
gray(seq(0,max(resmatrix_g),max(resmatrix_g)/(nrow(resmatrix_g)*ncol(resmatrix_g)))/max(resmatrix_g))


plot(c(0,(ncol(pheno)+1)),c(0.5,4.5), xlab="", ylab="", t="n")
for(a in 1:ncol(pheno)){
	rect((a-0.5),0.5,(a+0.5),1.5,col=gray(1-(0:(max(res)-min(res)+2))/(max(res)-min(res)+2))[(round(res[1,a]-min(res))+1)],border = "transparent")
	rect((a-0.5),1.5,(a+0.5),2.5,col=gray(1-(0:(max(res)-min(res)+2))/(max(res)-min(res)+2))[(round(res[2,a]-min(res))+1)],border = "transparent")
	rect((a-0.5),2.5,(a+0.5),3.5,col=gray(1-(0:(max(res)-min(res)+2))/(max(res)-min(res)+2))[(round(res[3,a]-min(res))+1)],border = "transparent")
	rect((a-0.5),3.5,(a+0.5),4.5,col=gray(1-(0:(max(res)-min(res)+2))/(max(res)-min(res)+2))[(round(res[4,a]-min(res))+1)],border = "transparent")
}



gray((0:(max(res)-min(res)+2))/(max(res)-min(res)+2))[(round(res[,a]-min(res))+1)]
col=rgb((0:(max(res)-min(res)+2))/(max(res)-min(res)+2),1-(0:(max(res)-min(res)+2))/(max(res)-min(res)+2),0.5)[(round(res[,a]-min(res))+1)]
