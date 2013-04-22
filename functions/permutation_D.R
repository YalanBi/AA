#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 15-04-2013
# first written: 15-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]


#direction selection
probesDir <- function(exp_data = rawexp){
  if(unique(exp_data[,"strand"]) == "sense"){
    direction_id <- which(exp_data[, "direction"] == "reverse")
  }
  if(unique(exp_data[,"strand"]) == "complement"){
    direction_id <- which(exp_data[, "direction"] == "forward")
  }
  return(direction_id)
}


#*************************************************************** load part ***************************************************************#
load(file="Data/fullModeMapping/expGenes.Rdata")


#************************************************************ make matrix part ************************************************************#
#make new expression files for permutation, for expressed genes, 1 tu 1 probe, at exon level
for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  genenames <- expGeneList[[chr]]
  
  newprobematrix <- NULL
  rownameList <- NULL
  colnameList <- NULL
  
  #filename="AT1G01010.txt"
  for(filename in genenames){
    rawexp <- read.table(paste0(location, filename), row.names=1, header=T)
    
    if(length(colnameList) == 0){
      colnameList <- c("bp", colnames(rawexp)[17:164])
      cat("get colnames!\n")
    }
    
    probes_dir <- probesDir(rawexp)
    exonID <- probes_dir[which(grepl("tu", rawexp[probes_dir, "tu"]))]
    
    #cat(filename)
    for(tu in unique(rawexp[exonID,"tu"])){
      rownameList <- c(rownameList, paste0(gsub(".txt", "", filename), "_", as.character(tu)))
      
      probes4tu <- exonID[which(rawexp[exonID,"tu"] == tu)]
      
      newprobematrix <- rbind(newprobematrix, c(round(mean(rawexp[probes4tu, "bp"])), colMeans(rawexp[probes4tu,17:164])))
      #cat(" ", as.character(tu))
    }
    #cat("\n")
  }
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s!\n\n")
  
  rownames(newprobematrix) <- rownameList
  colnames(newprobematrix) <- colnameList
  write.table(newprobematrix, file=paste0("Data/fullModeMapping/expGenes_chr", chr, "_1tu1probe.txt"), sep="\t")
}

#to load this file
#read.table("Data/fullModeMapping/expGenes_1gene1probe.txt", row.names=1, header=T)





#************************************************************ permutation part ************************************************************#
setwd("D:/Arabidopsis Arrays")
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]

#function for full model mapping
map.fast <- function(geno, pheno, menvironment){
  res <- NULL
  modelinfo <- summary(aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno) + as.factor(menvironment):as.numeric(geno)))
  res$qtl   <- unlist(lapply(modelinfo,"[",2,5),use.names=T)
  res$int   <- unlist(lapply(modelinfo,"[",3,5),use.names=T)
  res
}

#full model mapping with expressed genes
mapGenotypes <- function(geno, menvironment, permutation){
  st <- proc.time()
  cat("Start scanning for permutation", permutation, "\n")

  resQTL <- NULL
  resInt <- NULL
  
  for(m in 1:ncol(geno)){
    resList <- map.fast(geno[ ,m], pheno = expdata, menvironment)
    resQTL <- cbind(resQTL, -log10(resList$qtl))
    resInt <- cbind(resInt, -log10(resList$int))
  }
 
  et <- proc.time()
  cat("Done with full model mapping after:",(et-st)[3],"secs\n")
  return(c(max(resQTL), max(resInt)))
}

set.seed(1001)
#full model mapping, 5 chr
fileqtl <- "perms_qtls.txt"
fileint <- "perms_ints.txt"

cat("", file=fileqtl)
cat("", file=fileint)

chr <- 1
#the 1st col is bp!!!
expdata <- t(read.table(paste("Data/fullModeMapping/expGenes_chr", chr, "_1tu1probe.txt", sep=""),row.names=1, header=TRUE)[, 2:149])
expdata <- expdata[,1:100] # REMOVE BEFORE REAL !!!
  
for(permutation in 1:1000){
  newGeno <- geno[sample(1:148), ]
  scores <- mapGenotypes(geno = newGeno, menvironment, permutation)
  cat(paste0(permutation, "\t", scores[1],"\n"), file=fileqtl, append = TRUE)
  cat(paste0(permutation, "\t", scores[2],"\n"), file=fileint, append = TRUE)
}

