#P3-Generate images and matrix

#1---match genotypes in the order of v9-v172
#~~1) 2 for loops
setwd("X:/Data/Desktop/Arabidopsis Arrays")
geno <- read.table("genotypes.txt", colClasses="character")
ann <- read.table("annotation_sample.txt", colClasses="character")
geno_n <- NULL
ind_names_v <-  NULL
ind_names_r <-  NULL
for(x in 9:172){  #Line number of the annotation file
	for(y in 1: nrow(geno)){
		name <- paste("RIL", geno[y,1], sep="")
		if(ann[x,2] == name){
			ind_names_v <- c(ind_names_v, ann[x,1])
			ind_names_r <- c(ind_names_r, name)
			geno_n <- rbind(geno_n, as.character(geno[y,5:73]))
		}
	}
}
#~~~~check~~~~#
rownames(geno_n) <- ind_names_v
geno_n[1:10,1:10]
rownames(geno_n) <- ind_names_r
geno_n[1:10,1:10]

#~~2) 1 for loop + 1 which
geno_n <- NULL
ind_names_v <-  NULL
ind_names_r <-  NULL
geno_names <- paste("RIL", geno[,1], sep="")
for(x in 9:172){  #Line number of the annotation file
	if(ann[x,2] %in% geno_names){
		y <- which(geno_names %in% ann[x,2]) #"A %in% B" vs "B %in% A"#
		ind_names_v <- c(ind_names_v, ann[x,1])
		ind_names_r <- c(ind_names_r, geno_names[y])
		geno_n <- rbind(geno_n, as.character(geno[y,5:73]))
	}
}
write.table(geno_n, file="genotypes_n.txt", row.names = ind_names_r, col.names = FALSE)
#generate a new genotypes file!

#~~3) match
geno <- read.table("genotypes.txt", colClasses="character", row.names=paste("RIL", geno[,1], sep=""))
ann <- read.table("annotation_sample.txt", colClasses="character")
geno_n <- geno[match(ann[9:172,2], rownames(geno)),5:73]
write.table(geno_n, file="genotypes_n.txt", col.names = FALSE)




###2---do t-test for genetic effects###
geno__n <- read.table("genotypes_n.txt", row.names=1)
for(filename in dir("Data/gene_data")[grepl(".txt",dir("Data/gene_data"))]){
	cat("Loading", filename,"\n")
	Pheno <- read.table(paste("Data/gene_data/",filename,sep=""))
	pheno <- t(Pheno[,25:188])
	resmatrix <- NULL
	png(file = paste("Data/QTL Profile~t.test/",gsub(".txt",".png",filename),sep=""), bg="white")
	plot(c(0,69), c(0,20), t="n")
	for(p in 1:ncol(pheno)){
		genet_eff <- NULL
		for(x in 1:ncol(geno__n)){
			g1 <- which(geno__n[,x] == "A")
			g2 <- which(geno__n[,x] == "B")
			p_value <- t.test(pheno[g1,p],pheno[g2,p])$p.value
			genet_eff <- c(genet_eff, -log(p_value))
		}
		points(genet_eff, t='l', col=p)
	}
	dev.off()
	return()
}



#2---do lm for genetic effect
setwd("X:/Data/Desktop/Arabidopsis Arrays")

geno <- read.table("genotypes_n.txt", row.names=1)
menvironment <- read.table("annotation_sample.txt", colClasses="character")[9:172,3]

#~~1) anova
for(filename in dir("Data/gene_data")[grepl(".txt",dir("Data/gene_data"))]){
	cat("Loading", filename,"\n")
	Pheno <- read.table(paste("Data/gene_data/",filename,sep=""))
	pheno <- t(Pheno[,25:188])
	resmatrix <- NULL
	for(p in 1:ncol(pheno)){
		genet_eff <- NULL
		for(x in 1:ncol(geno)){
			mylm <- anova(lm(pheno[,p] ~ as.factor(menvironment) + as.numeric(geno[,x])))
			genet_eff <- c(genet_eff, -log(mylm[[5]][2]))
		}
		resmatrix <- rbind(resmatrix, genet_eff)
	}
	write.table(resmatrix, file=paste("Data/gene_eff~lm/", gsub(".txt","_G.txt",filename), sep=""), row.names = colnames(pheno), col.names = paste("g", 1:69, sep=""))
	return()
}

#~~2) aov
map.fast <- function(x, geno, pheno, menvironment){
	res <- NULL
	models  <- aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno[,x]))
	modelinfo <- summary(models)
	res$env <- unlist(lapply(modelinfo,"[",1,5),use.names=T)
	res
}
for(filename in dir("Data/gene_data")[grepl(".txt",dir("Data/gene_data"))]){
	cat("Loading", filename,"\n")
	Pheno <- read.table(paste("Data/gene_data/",filename,sep=""))
	pheno <- t(Pheno[,25:188])
	resmatrix_g <- NULL
	for(x in 1:ncol(geno)){
		resmatrix_g <- cbind(resmatrix_g, -log(map.fast(x, geno, pheno, menvironment)$qtl))
	}
	write.table(resmatrix_g, file=paste("Data/gene_eff~lm/", gsub(".txt","_G.txt",filename), sep=""), row.names = colnames(pheno), col.names = paste("g", 1:69, sep=""))
	#return()
}

#3---images of lm~gene_eff after cut-off at 4
for(filename in dir("Data/gene_eff~lm")[grepl(".txt",dir("Data/gene_eff~lm"))]){
	cat("Loading", filename,"\n")
	resmatrix_g <- read.table(paste("Data/gene_eff~lm/",filename,sep=""))
	png(file = paste("Data/QTL Profile~lm/",gsub("_G.txt","~lm.png",filename),sep=""), bg="white")
	plot(c(0,69), c(0,20), t="n")
	for(g in 1:nrow(resmatrix_g)){
		points(seq(1,69),resmatrix_g[g,], t='l', col=g)
	}
	dev.off()
	#return()
}

#4---summerize the probes into groups
for(filename in dir("Data/gene_eff~lm")[grepl(".txt",dir("Data/gene_eff~lm"))]){
  gene_eff <- read.table(paste("Data/gene_eff~lm/",filename,sep=""))
  for(g in 1:ncol(gene_eff)){
    g_cut <- as.numeric(gene_eff[,g] >= 4)
    if(sum(g_cut) >= 4){
      group <- which(g_cut == 1)
      cat(gsub("G.txt",paste("m", g, ":", sep=""),filename), group, "\n", sep=" ", file="group_ind~cut-off=4.txt", append=TRUE)
    }
  }
  cat("\n", file="group_ind~cut-off=4.txt", append=TRUE)
  #return()
}


for(filename in dir("Data/gene_eff~lm")[grepl(".txt",dir("Data/gene_eff~lm"))]){
	gene_eff <- read.table(paste("Data/gene_eff~lm/",filename,sep=""))
	for(g in 1:ncol(gene_eff)){
		g_cut <- as.numeric(gene_eff[,g] >= 4)
		if(sum(g_cut) > 0){
			group <- which(g_cut == 1)
			cat(gsub("G.txt",paste("m", g, ":", sep=""),filename), group, "\n", sep=" ", file="group_ind~cut-off=4_2.txt", append=TRUE)
		}
	}
	cat("\n", file="group_ind~cut-off=4_2.txt", append=TRUE)
	#return()
}










#images of correlation
image_cor <- function(){
	png(file = paste("image_lm_cor/",gsub(".txt","_cor.png",filename),sep=""), bg="white")
	image(cor(resmatrix))
	dev.off()
}
#images of heatmap
image_hm <- function(){
	png(file = paste("image_lm_heatmap/",gsub(".txt","_hm.png",filename),sep=""), bg="white", width = 960, height = 960,)
	heatmap(cor(resmatrix),scale="none")
	dev.off()
}


for(filename in dir("gene_data")[grepl(".txt",dir("gene_data"))]){
	cat("Loading", filename,"\n")
	Pheno <- read.table(paste("gene_data/",filename,sep=""))
	pheno <- t(Pheno[,25:188])
	resmatrix1 <- NULL
	#~~images of lm~gene_eff after cut-off at 4
	#~~png(file = paste("gene_data/image_lm_black/",gsub(".txt","_b.png",filename),sep=""), bg="white")
	#~~plot(c(0,69), c(0,(ncol(pheno)+1)), t="n")
	for(p in 1:ncol(pheno)){
		genet_eff <- NULL
		for(x in 1:ncol(geno)){
			mylm <- anova(lm(pheno[,p] ~ as.factor(menvironment) + as.numeric(geno[,x])))#############aov <- Function with a for loop##############
			genet_eff <- c(genet_eff, -log(mylm[[5]][2]))
		}
		#~~points(as.numeric(genet_eff >= 4)*p, t='p', col="black", pch=20,cex=3)
		#cat("Pheno", p, ": ", genet_eff, "\n")
		resmatrix1 <- rbind(resmatrix, genet_eff)
	}
	#~~dev.off()
	#write.table(resmatrix, file=paste("genet_eff/", gsub(".png",".txt",filename), sep=""), row.names = colnames(pheno), col.names = paste("g", 1:69, sep=""))
	for(g in 1: ncol(geff_cut)){
		if(sum(geff_cut[,g]) >= 4){
			probes <- NULL
			if(geff_cut[,g]){
				
			}
			c(filename, "_", a, ": ", )
		}
	}
	row.names(resmatrix) = colnames(pheno)
	if(nrow(resmatrix) != 1){
		#image_cor()
		#image_hm()
	}
	return()
}








#~~~test for lm~~~
setwd("X:/Data/Desktop/Arabidopsis Arrays/Original Data")
#setwd("/Users/SiBYL/Desktop/Arabidopsis Arrays/Original Data")
geno <- read.table("genotypes.txt", colClasses="character")
ann <- read.table("annotation_sample.txt", colClasses="character")
geno__n <- read.table("genotypes_n.txt", row.names=1)

geno_n <- NULL
ind_names_v <-  NULL
ind_names_r <-  NULL

menvironment <- ann[9:172,3]



#images of correlation
image_cor <- function(){
	png(file = paste("image_lm_cor/",gsub(".txt","_cor.png",filename),sep=""), bg="white")
	image(cor(resmatrix))
	dev.off()
}
#images of heatmap
image_hm <- function(){
	png(file = paste("image_lm_heatmap/",gsub(".txt","_hm.png",filename),sep=""), bg="white", width = 960, height = 960,)
	heatmap(cor(resmatrix),scale="none")
	dev.off()
}
###~~~daiding~~~###
genet_eff <- function(p){
	g_eff <- NULL
	for(x in 1:ncol(geno__n)){
		mylm <- anova(lm(pheno[,p] ~ as.factor(menvironment) + as.numeric(geno__n[,x])))#############aov <- Function with a for loop##############
		g_eff <- c(g_eff, -log(mylm[[5]][2]))
	}
	return(g_eff)
}
genet_eff(1)


for(filename in dir("gene_data")[grepl(".txt",dir("gene_data"))]){
	cat("Loading", filename,"\n")
	Pheno <- read.table(paste("gene_data/",filename,sep=""))
	pheno <- t(Pheno[,25:188])
	resmatrix <- NULL
	#~~images of lm~gene_eff after cut-off at 4
	#~~png(file = paste("gene_data/image_lm_black/",gsub(".txt","_b.png",filename),sep=""), bg="white")
	#~~plot(c(0,69), c(0,(ncol(pheno)+1)), t="n")
	for(p in 1:ncol(pheno)){
		genet_eff <- NULL
		for(x in 1:ncol(geno__n)){
			mylm <- anova(lm(pheno[,p] ~ as.factor(menvironment) + as.numeric(geno__n[,x])))#############aov <- Function with a for loop##############
			genet_eff <- c(genet_eff, -log(mylm[[5]][2]))
		}
		#~~points(as.numeric(genet_eff >= 4)*p, t='p', col="black", pch=20,cex=3)
		#cat("Pheno", p, ": ", genet_eff, "\n")
		resmatrix <- rbind(resmatrix, genet_eff)
	}
	#~~dev.off()
	#write.table(resmatrix, file=paste("genet_eff/", gsub(".png",".txt",filename), sep=""), row.names = colnames(pheno), col.names = paste("g", 1:69, sep=""))
	for(g in 1: ncol(geff_cut)){
		if(sum(geff_cut[,g]) >= 4){
			probes <- NULL
			if(geff_cut[,g]){
				
			}
			c(filename, "_", a, ": ", )
		}
	}
	row.names(resmatrix) = colnames(pheno)
	if(nrow(resmatrix) != 1){
		#image_cor()
		#image_hm()
	}
	return()
}

#make a matrix <- 0,1;39*69;in a col, sum>=4 -> a group
		geff_cut <- rbind(geff_cut, as.numeric(genet_eff >= 4))
#make a file <- groups
#cor(expression data) 164*164;cut-off = mean;a group
#cor(profile: gene_eff data) 69*69;cut-off = mean;a group


#Go through all result matrices and create groups
