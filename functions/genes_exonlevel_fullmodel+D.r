#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 23-04-2013
# first written: 19-04-2013
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
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
    
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





#*************************************************** full model mapping with directions ***************************************************#
setwd("D:/Arabidopsis Arrays")
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]

#function for full model mapping
map.fast <- function(geno, pheno, menvironment){
  res <- NULL
  models    <- aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno) + as.factor(menvironment):as.numeric(geno))
  modelinfo <- summary(models)
  res$env   <- unlist(lapply(modelinfo,"[",1,5),use.names=T)
  res$qtl   <- unlist(lapply(modelinfo,"[",2,5),use.names=T)
  res$int   <- unlist(lapply(modelinfo,"[",3,5),use.names=T)
  res$eff   <- unlist(models$coefficients[2,])
  res
}

#full model mapping with expressed genes
mapGenotypes <- function(chr=1, geno, menvironment){
  st <- proc.time()
  cat("Loading chr", chr, "\n")
  
  resEnv <- NULL
  resQTL <- NULL
  resInt <- NULL
  
  for(m in 1:ncol(geno)){
    resList <- map.fast(geno[ ,m], pheno = expdata, menvironment)
    resEnv <- cbind(resEnv, -log10(resList$env))
    resQTL <- cbind(resQTL, -log10(resList$qtl) * sign(resList$eff))
    resInt <- cbind(resInt, -log10(resList$int))
  }
  
  rownames(resEnv) <- colnames(expdata)
  colnames(resEnv) = colnames(geno)
  write.table(resEnv, file=paste("Data/fullModeMapping/expGenes_chr", chr, "_FMD_Env.txt", sep=""), sep="\t")
  
  rownames(resQTL) <- colnames(expdata)
  colnames(resQTL) = colnames(geno)
  write.table(resQTL, file=paste("Data/fullModeMapping/expGenes_chr", chr, "_FMD_QTL.txt", sep=""), sep="\t")
  
  rownames(resInt) <- colnames(expdata)
  colnames(resInt) = colnames(geno)
  write.table(resInt, file=paste("Data/fullModeMapping/expGenes_chr", chr, "_FMD_Int.txt", sep=""), sep="\t")
  
  et <- proc.time()
  cat("Done with full model mapping after:",(et-st)[3],"secs\n")
}

#full model mapping, 5 chr
for(chr in 1:5){
  #the 1st col is bp!!!
  expdata <- t(read.table(paste("Data/fullModeMapping/expGenes_chr", chr, "_1tu1probe.txt", sep=""),row.names=1, header=TRUE)[, 2:149])
  mapGenotypes(chr = chr, geno, menvironment)
}
