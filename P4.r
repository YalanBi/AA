#P4-summery for chromosome
setwd("X:/Data/Desktop/Arabidopsis Arrays")
ann_m <- t(read.table("annotation_marker.txt"))

for(filename in dir("Data/gene_eff~lm")[grepl(".txt",dir("Data/gene_eff~lm"))]){
	res <- group(filename, lodThreshold = 4)
	allprobes <- 1:res$nprobes
	qtlprobes <- unique(unlist(res$ind))
	expprobes <- which(res$means > 2)
	goodprobes <- unique(c(qtlprobes,expprobes))
	badprobes <- which(!allprobes %in% goodprobes)
	cat(filename, badprobes,"\n")
}


st  <- proc.time()
ratiomatrix <- NULL
for(filename in dir("Data/gene_eff~lm")[grepl(".txt",dir("Data/gene_eff~lm"))]){
  cat("loading",filename,"\n")
  #cat(gsub("_G.txt","",filename), "\n", sep="", file="Data/ratio_ID.txt", append=T)
  ratiomatrix <- rbind(ratiomatrix,group(filename, lodThreshold = 4)$ratio)
  #for(chr in 1:5){
    #cat(c(group(filename, lodThreshold = 4)$ratio[chr],as.character(group(filename, lodThreshold = 4)$ind[[chr]])), "\n", file="Data/ratio_ID.txt", append=T)
  #}
  #cat("\n\n", sep="", file="Data/ratio_ID.txt", append=T)
  return()
}
rownames(ratiomatrix) <- gsub("_G.txt","",dir("Data/gene_eff~lm")[grepl(".txt",dir("Data/gene_eff~lm"))])
colnames(ratiomatrix) <- paste("chr",1:5,sep="")
write.table(ratiomatrix, file="Data/ratio.txt")
cat("Done in:",(proc.time()-st)[3],"seconds\n")



st  <- proc.time()
cat(c("\t", "chr","maxratio","probes_ind"), "\n", file="Data/ratio_max.txt", append=T)
for(filename in dir("Data/gene_eff~lm")[grepl(".txt",dir("Data/gene_eff~lm"))]){
  cat("loading",filename,"\n")
  if(group(filename, lodThreshold = 4)$ratio[which.max(group(filename, lodThreshold = 4)$ratio)] != 0){
  cat(c(gsub("_G.txt","",filename), which.max(group(filename, lodThreshold = 4)$ratio), group(filename, lodThreshold = 4)$ratio[which.max(group(filename, lodThreshold = 4)$ratio)],as.character(group(filename, lodThreshold = 4)$ind[[which.max(group(filename, lodThreshold = 4)$ratio)]])), "\n", file="Data/ratio_max.txt", append=T)
  }
  return()
}
cat("Done in:",(proc.time()-st)[3],"seconds\n")




group(filename = "Data/gene_data/AT1G01115.txt", lodThreshold = 4)














getChromo <- function(annotation, chr=1){
	which(as.numeric(t(annotation)[,2])==chr)
}



setwd("X:/Data/Desktop/Arabidopsis Arrays")
annotation <- read.table("annotation_marker.txt",sep="\t")






image(y=1:ncol(resmatrix_g),x=1:nrow(resmatrix_g), z=resmatrix_g)
grid(nrow(resmatrix_g),ncol(resmatrix_g),col="white",lty=1)
box()