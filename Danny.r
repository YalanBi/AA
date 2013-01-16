setwd("X:/Data/Desktop/Arabidopsis Arrays")
geno <- read.table("Data/genotypes_n.txt", row.names=1)
ann <- read.table("annotation_sample.txt", colClasses="character")
menvironment <- ann[9:172,3]
annotation <- read.table("annotation_marker.txt")
map.fast <- function(x, geno, pheno, menvironment){
	res <- NULL
	models  <- aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno[,x]))
	modelinfo <- summary(models)
	res$env <- unlist(lapply(modelinfo,"[",1,5),use.names=T)
	res$qtl <- unlist(lapply(modelinfo,"[",2,5),use.names=T)
	res
}
getChromo <- function(annotation, chr=1){
	which(as.numeric(t(annotation)[,2])==chr)
}

getRatioOfProbesWithQTL <- function(filename = "Data/gene_data/AT1G01010.txt", lodThreshold = 4){
	expressiondata <- read.table(filename)
	resmatrix_g <- NULL
	for(x in 1:ncol(geno)){
		resmatrix_g <- cbind(resmatrix_g, -log(map.fast(x, geno, t(expressiondata[,25:188]), menvironment)$qtl))
	}
	v <- NULL
	idsmatrix <- vector("list",5)
	for(chr in 1:5){
		s <- 0
		ids <- NULL
		for(x in 1:nrow(resmatrix_g)){
			if(any(resmatrix_g[x ,getChromo(annotation, chr)] > lodThreshold)){
				s <- s + 1
				ids  <- c(ids, x)
			}
		}
		idsmatrix[[chr]] <- ids
		v <- c(v,s)
	}
	return(list(v/nrow(resmatrix_g), idsmatrix))
}
ratiomatrix <- NULL
for(filename in dir("Data/gene_data")[grepl(".txt",dir("Data/gene_data"))]){
  cat("loading",filename,"\n")
  ratiomatrix <- rbind(ratiomatrix,c(gsub(".txt","",filename),getRatioOfProbesWithQTL(paste("Data/gene_data/",filename, sep=""), lodThreshold = 4)[[1]]))
  #return()
}
colnames(ratiomatrix) <- c("genenames",paste("chr",1:5,sep=""))
write.table(ratiomatrix, file="Data/ratio.txt", row.names=FALSE)
read.table("Data/ratio.txt")