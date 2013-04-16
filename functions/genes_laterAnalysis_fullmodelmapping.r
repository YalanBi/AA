#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 15-04-2013
# first written: 15-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")


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


#*************************************************************** test part ***************************************************************#
#find expressed genes, median of all exon RILs of this gene higher than threshold=5
findExpGene <- function(chr, threshold = 5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  expGene <- NULL
  genenames <- dir(location)[which(grepl(".txt", dir(location)) & !grepl("_QTL", dir(location)))]
  
  #filename "AT1G01010.txt"
  for(filename in genenames){
    rawexp <- read.table(paste0(location, filename), row.names=1, header=T)
    
    probes_dir <- probesDir(rawexp)
    exonID <- probes_dir[which(grepl("tu", rawexp[probes_dir, "tu"]))]
    
    if(length(exonID) > 0){
      cat(filename, "median is", median(unlist(rawexp[exonID, 17:164])))
      if(median(unlist(rawexp[exonID, 17:164])) >= threshold){
        expGene <- c(expGene, filename)
        cat(" TRUE!")
      }
      cat("\n")
    }
  }
  et <- proc.time()[3]
  cat("finding expressed genes on chr", chr, "finished in", et-st, "s\n")
  
  return(expGene)
}
#expGene <- findExpGene(chr = 1, threshold = 5)


#*************************************************************** save part ***************************************************************#
#if necessary, save the results into a file
expGeneList <- vector("list", 5)
for(chr in 1:5){
  expGeneList[[chr]] <- findExpGene(chr, threshold = 5)
}
save(expGeneList,  file="Data/fullModeMapping/expGenes.Rdata")


#************************************************************** mapping part **************************************************************#
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]

#function for full model mapping
map.fast <- function(geno, pheno, menvironment){
  res <- NULL
  modelinfo <- summary(aov(as.matrix(pheno) ~ as.factor(menvironment) + as.numeric(geno) + as.factor(menvironment):as.numeric(geno)))
  res$env   <- unlist(lapply(modelinfo,"[",1,5),use.names=T)
  res$qtl   <- unlist(lapply(modelinfo,"[",2,5),use.names=T)
  res$int   <- unlist(lapply(modelinfo,"[",3,5),use.names=T)
  res
}

#full model mapping with expressed genes
mapGenotypes <- function(filename, geno, menvironment, chr=1){
  #filename=AT1G01010.txt
  ###only check Int file!!###
  if(!file.exists(paste("Data/fullModeMapping/chr", chr, "_norm_hf_cor/", gsub(".txt", "_FM_Int.txt", filename), sep=""))){
    st <- proc.time()
    cat("Loading", filename, "\n")
    
    rawexp <- read.table(paste("Data/chr", chr, "_norm_hf_cor/", filename, sep=""),row.names=1, header=TRUE)
    newexp <- t(rawexp[,17:ncol(rawexp)])#####col numbers are changed!!!#####
    resList <- NULL
    resEnv <- NULL
    resQTL <- NULL
    resInt <- NULL
    
    for(m in 1:ncol(geno)){
      resList <- map.fast(geno[ ,m], pheno = newexp, menvironment)
      resEnv <- cbind(resEnv, -log10(resList$env))
      resQTL <- cbind(resQTL, -log10(resList$qtl))
      resInt <- cbind(resInt, -log10(resList$int))
    }
    
    rownames(resEnv) <- colnames(newexp)
    colnames(resEnv) = colnames(geno)
    write.table(resEnv, file=paste("Data/fullModeMapping/chr", chr, "_norm_hf_cor/", gsub(".txt","_FM_Env.txt",filename), sep=""))
    
    rownames(resQTL) <- colnames(newexp)
    colnames(resQTL) = colnames(geno)
    write.table(resQTL, file=paste("Data/fullModeMapping/chr", chr, "_norm_hf_cor/", gsub(".txt","_FM_QTL.txt",filename), sep=""))
    
    rownames(resInt) <- colnames(newexp)
    colnames(resInt) = colnames(geno)
    write.table(resInt, file=paste("Data/fullModeMapping/chr", chr, "_norm_hf_cor/", gsub(".txt","_FM_Int.txt",filename), sep=""))
    
    et <- proc.time()
    cat("Done with full model mapping after:",(et-st)[3],"secs\n")
  } else{
    cat("Skipping", filename,", because it exists\n")
  }
}

#full model mapping, 5 chr
for(chr in 1:5){
  expGene <- findExpGene(chr, threshold = 5)
  
  for(filename in expGene){
    mapGenotypes(filename, geno, menvironment, chr)
  }
}
