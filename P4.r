#P4-summery for chromosome
setwd("X:/Data/Desktop/Arabidopsis Arrays")
geno <- read.table("Data/genotypes_n.txt", row.names=1)
ann <- read.table("annotation_sample.txt", colClasses="character")
menvironment <- ann[9:172,3]
ann_m <- t(read.table("annotation_marker.txt"))

group_chr <- function(ann_m, chr){
	which(ann_m[,2] == chr)
}
group <- function(filename, lodThreshold = 4){
	gene_eff <- read.table(paste("Data/gene_eff~lm/",filename,sep=""))
	sum_num <- NULL
	indmatrix <- NULL
	res <- NULL
	for(chr in 1:5){
		s <- 0
		ind <- NULL
		for(p in 1:nrow(gene_eff)){
			if(any(gene_eff[p,group_chr(ann_m, chr)] >= lodThreshold)){
				s <- s+1
				ind <- c(ind, p)
			}
		}
		 sum_num <- c(sum_num, s)
		 indmatrix[[chr]] <- ind
	}
	res$ratio <- sum_num/nrow(gene_eff)
	res$ind <- indmatrix
	cat("res$ratio","\n",file="resuls.txt",append=T)
	cat(res$ratio,"\n",file="resuls.txt",append=T)
	cat("res$ind","\n",file="resuls.txt",append=T)
	cat(res$ind,"\n",file="resuls.txt",append=T) ###?????###
	return(res)
}

for(filename in dir("Data/gene_eff~lm")[grepl(".txt",dir("Data/gene_eff~lm"))]){
	group(filename, lodThreshold = 4)
}
$ratio
[1] 0.20512821 0.07692308 0.12820513 0.12820513 0.30769231

$ind
$ind[[1]]
[1]  1  4  7 10 16 17 24 31

$ind[[2]]
[1] 12 30 33

$ind[[3]]
[1]  6  9 26 28 37

$ind[[4]]
[1]  9 14 15 28 35

$ind[[5]]
 [1]  1  2 12 14 15 19 24 26 27 29 37 39




















getChromo <- function(annotation, chr=1){
	which(as.numeric(t(annotation)[,2])==chr)
}



setwd("X:/Data/Desktop/Arabidopsis Arrays")
annotation <- read.table("annotation_marker.txt",sep="\t")






image(y=1:ncol(resmatrix_g),x=1:nrow(resmatrix_g), z=resmatrix_g)
grid(nrow(resmatrix_g),ncol(resmatrix_g),col="white",lty=1)
box()