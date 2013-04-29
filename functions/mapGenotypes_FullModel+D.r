#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 29-04-2013
# first written: 29-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]

#*************************************************** full model mapping with directions ***************************************************#
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

mapAndSaveQTL <- function(geno, pheno, menvironment, location, filename, fstem = "_1tu1probe.txt"){
  st <- proc.time()
  resList <- apply(geno, 2, function(x, pheno, env){
    result <- map.fast(x, pheno, menvironment)
    return(list(-log10(result$env), sign(result$eff) * -log10(result$qtl), -log10(result$int)))
  }, pheno=rawexp, env=menvironment)
  
  qtls <- lapply(resList,"[[",2)
  qtlm <- matrix(unlist(qtls), length(qtls[[1]]), length(qtls))
  rownames(qtlm) <- colnames(rawexp)
  colnames(qtlm) <- colnames(geno)
  write.table(qtlm, file=paste0(location, gsub(fstem, "_FMD_QTL_Danny.txt", filename)), sep="\t") #NOTICE: change name!!!
  
  envs <- lapply(resList,"[[",1)
  envm <- matrix(unlist(envs), length(envs[[1]]), length(envs))
  rownames(envm) <- colnames(rawexp)
  colnames(envm) <- colnames(geno)
  write.table(envm, file=paste0(location, gsub(fstem, "_FMD_Env_Danny.txt", filename)), sep="\t") #NOTICE: change name!!!
  
  ints <- lapply(resList,"[[",3)
  intm <- matrix(unlist(ints), length(ints[[1]]), length(ints))
  rownames(intm) <- colnames(rawexp)
  colnames(intm) <- colnames(geno)
  write.table(intm, file=paste0(location, gsub(fstem, "_FMD_Int_Danny.txt", filename)), sep="\t") #NOTICE: change name!!!

  et <- proc.time()
  cat("Done with QTL mapping after:",(et-st)[3],"secs\n")
}




#************************************************************** mapping part **************************************************************#
#exon level: 1tu1probe, 5 files
location_FM <- "Data/fullModeMapping/"
for(chr in 1:5){
  filename <- paste0("expGenes_chr", chr, "_1tu1probe.txt")
  rawexp <- t(read.table(paste(location_FM, filename, sep=""),row.names=1, header=TRUE)[, 2:149])
  mapAndSaveQTL(geno, pheno = rawexp, menvironment, location = location_FM, filename, fstem = "_1tu1probe.txt")
}




#probe level: need to be updated when used...
for(filename in dir(paste("Data/chr", chr, "_norm_hf_cor/", sep=""))[grepl(".txt",dir(paste("Data/chr", chr, "_norm_hf_cor/", sep=""))) & !grepl("_QTL",dir(paste("Data/chr", chr, "_norm_hf_cor/", sep="")))]){
  if(!file.exists(paste("Data/chr", chr, "_norm_hf_cor_fullModel/", gsub(".txt","_FM_QTL.txt",filename), sep=""))){ ###only check QTL file!!###
    st <- proc.time()
    cat("Loading", filename,"\n")
    Pheno <- read.table(paste("Data/chr", chr, "_norm_hf_cor/", filename, sep=""),row.names=1, header=TRUE)
    if(nrow(Pheno) < 4){
      cat("Skipping", filename," because it has to few probes\n")
    }else{
      pheno <- t(Pheno[,17:ncol(Pheno)])#####col numbers are changed!!!#####
      mapAndSaveQTL(geno, pheno, menvironment, chr, filename)
    }
  }else{
    cat("Skipping", filename," because it exists\n")
  }
}
