#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 26-03-2013
# first written: 26-03-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#
###########
# library #
###########

setwd("D:/Arabidopsis Arrays/Data/newTU")

###################################### For chr1&5(_t) loding into R ######################################
getGene <- function(chr = 1, what="QTL", startG = "AT1G01010", untill = "AT1G21835"){
  x <- 0
  matrix <- NULL
  mnames <- NULL
  fp <- file(paste0("genesByTU_chr", chr, "_norm_hf_cor_", what, "_t.txt"))
  open(fp)
  line <- readLines(fp, 1)
  tnames <- strsplit(line,"\t")[[1]][-1]
  line <- readLines(fp, 1)
  x <- x+1
  
  s <- min(which(grepl(startG, tnames)))
  e <- max(which(grepl(untill, tnames)))
  
  while(length(line) > 0 && line != ""){
    data <- strsplit(line,"\t")[[1]]
    mnames <- c(mnames, data[1])
    cat("Done",x,"\n")
    
    data <- round(as.numeric(data[-1]), d=2)[s:e]
    matrix <- cbind(matrix, data)
    line <- readLines(fp, 1)
    
    x <- x+1
  }
  colnames(matrix) <- mnames
  rownames(matrix) <- tnames[s:e]
  
  close(fp)
  invisible(matrix)
}

######### split into 2 parts #########
#chr1 - part1
exp_data <- read.table("genesByTU_chr1_norm_hf_cor.txt", row.names=1, header=T)
qtls <- NULL
ints <- NULL
envs <- NULL
qtls <- getGene(chr = 1, "QTL", startG = "AT1G01010", untill = "AT1G38185")
ints <- getGene(chr = 1, "Int", startG = "AT1G01010", untill = "AT1G38185")
envs <- getGene(chr = 1, "Env", startG = "AT1G01010", untill = "AT1G38185")
#chr1 - part2
exp_data <- read.table("genesByTU_chr1_norm_hf_cor.txt", row.names=1, header=T)
qtls <- NULL
ints <- NULL
qtls <- getGene(chr = 1, "QTL", startG = "AT1G38460", untill = "AT1G81020")
ints <- getGene(chr = 1, "Int", startG = "AT1G38460", untill = "AT1G81020")
envs <- getGene(chr = 1, "Env", startG = "AT1G38460", untill = "AT1G81020")


#chr5 - part1
exp_data <- read.table("genesByTU_chr5_norm_hf_cor.txt", row.names=1, header=T)
qtls <- NULL
ints <- NULL
qtls <- getGene(chr = 5, "QTL", startG = "AT5G01010", untill = "AT5G36960")
ints <- getGene(chr = 5, "Int", startG = "AT5G01010", untill = "AT5G36960")
envs <- getGene(chr = 5, "Env", startG = "AT5G01010", untill = "AT5G36960")

#chr5 - part2
exp_data <- read.table("genesByTU_chr5_norm_hf_cor.txt", row.names=1, header=T)
qtls <- NULL
ints <- NULL
qtls <- getGene(chr = 5, "QTL", startG = "AT5G36970", untill = "AT5G67630")
ints <- getGene(chr = 5, "Int", startG = "AT5G36970", untill = "AT5G67630")
envs <- getGene(chr = 5, "Env", startG = "AT5G36970", untill = "AT5G67630")




###################################### For chr2,3,4 loding into R ######################################
chr = 4

exp_data <- read.table(paste0("genesByTU_chr", chr, "_norm_hf_cor.txt"), row.names=1, header=T)
qtls <- read.table(paste0("genesByTU_chr", chr, "_norm_hf_cor_QTL.txt"), row.names=1, header=T)
ints <- read.table(paste0("genesByTU_chr", chr, "_norm_hf_cor_Int.txt"), row.names=1, header=T)
envs <- read.table(paste0("genesByTU_chr", chr, "_norm_hf_cor_Env.txt"), row.names=1, header=T)





createGeneObject <- function(name = "AT1G01010", exp_data, qtls, ints, chr = 1){
  gene_id <- which(grepl(name, rownames(exp_data)))
  ntus <- length(gene_id)
  
  res <- NULL
  res$Exp <- exp_data[gene_id,]
  for(t in 1:ntus){
    gene_tu <- rownames(exp_data)[gene_id][t]
    tuID <- gsub(paste0(name, "_"), "", gene_tu)
    res[[tuID]] <- vector("list")
    res[[tuID]]$QTL <- qtls[gene_tu,]
    res[[tuID]]$Int <- ints[gene_tu,]
    res[[tuID]]$Env <- envs[gene_tu,]
  }
  save(res, file=paste0("geneObject_chr", chr, "/", name, ".Rdata"))
}



for(name in unlist(unique(lapply(strsplit(rownames(qtls), "_"), "[[", 1)))){
  createGeneObject(name, exp_data, qtls, ints, chr = 4)#change chr!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  cat(name, "done\n")
}


res
length(res)
length(res$tu1)
rownames(res$Exp)



classifySingleGene <- function(geneObject){

}









##########################################################################################################
#find the middle gene of chr1
genenames <- unlist(unique(lapply(strsplit(rownames(exp_data), "_"), "[[", 1)))
length(genenames) #7242 unique genes in chr1
genenames[length(genenames)/2] #"AT1G38185"
which(grepl(genenames[length(genenames)/2], rownames(exp_data))) #15505
#max(which(grepl(genenames[length(genenames)/2], rownames(exp_data)))) #if have multiple tus
rownames(exp_data)[which(grepl(genenames[length(genenames)/2], rownames(exp_data)))] #"AT1G38185_tu1"
# middle + 1
genenames[length(genenames)/2 + 1] #"AT1G38460"
rownames(exp_data)[which(grepl(genenames[length(genenames)/2 + 1], rownames(exp_data)))] #"AT1G38460_tu1"
# last gene
genenames[length(genenames)] #"AT1G81020"
##########################################################################################################
#find the middle gene of chr5
genenames <- unlist(unique(lapply(strsplit(rownames(exp_data), "_"), "[[", 1)))
length(genenames) #6518 unique genes in chr1
genenames[length(genenames)/2] #"AT5G36960"
which(grepl(genenames[length(genenames)/2], rownames(exp_data))) #13539 13540
max(which(grepl(genenames[length(genenames)/2], rownames(exp_data)))) #13540     #if have multiple tus
rownames(exp_data)[max(which(grepl(genenames[length(genenames)/2], rownames(exp_data))))] #"AT5G36960_tu2"
# middle + 1
genenames[length(genenames)/2 + 1] #"AT5G36970"
rownames(exp_data)[which(grepl(genenames[length(genenames)/2 + 1], rownames(exp_data)))] #"AT5G36970_tu2"
# last gene
genenames[length(genenames)] #"AT5G67630"
##########################################################################################################