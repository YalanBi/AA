at1 <- read.table("Data/Chr1/AT1G01010.txt",header=TRUE,row.names=1)
ori <- read.table("Data/Old/raw data/exp_ann30K.txt",header=TRUE,row.names=1,nrows=500)
ann <- read.table("Data/Old/raw data/annotation_sample.txt")


for(x in as.character(ann[,2])){
  cat(x,ori[50,x],at1[1,x],"\n",sep="   ")
}


geno_ori <- read.table("Data/Old/raw data/genotypes.txt")
namez <- NULL
for(y in 1: nrow(geno_ori)){
  namez <- c(namez, paste("RIL", geno_ori[y,1], sep=""))
}
geno_ori <- geno_ori[,-c(1:4)]

rownames(geno_ori) <- namez


mar_ann <- read.table("Data/Old/raw data/annotation_marker.txt")

colnames(geno_ori) <- as.character(unlist(mar_ann[1,]))

geno_new <- read.table("refined map/genotypes.txt")

colnames(ori) <- c(colnames(ori)[1:16],as.character(ann[,2]))

geno_ord <- geno_new[,colnames(geno_ori)]

for(x in as.character(ann[,2])){
  cat(x,geno_ori[x,1],geno_ord[x,1],"\n",sep="   ")
}

corsss <-  NULL
for(x in 1:ncol(geno)){
  corsss <- c(corsss, cor(as.numeric(geno[,x]),geno_new[,326]))
}