#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 18-04-2013
# first written: 18-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
#load SNP probes' number, snp is a list of 475009 integers
snp <- read.table("snp ids.txt")[,1]


#*************************************************************** remove part ***************************************************************#
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts now!\n")
  
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  location_FM <- paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/")
  location_snpCor <- paste0("Data/snpCorrection/chr", chr, "_norm_hf_cor/")
  
  #filename = AT1G01010.txt
  for(filename in dir(location)[which(grepl(".txt", dir(location)) & !grepl("_QTL", dir(location)) & !grepl("_SNP", dir(location)))]){
    rawexp <- read.table(paste0(location, filename), row.names=1, header=T)
    
    if(any(rownames(rawexp) %in% snp)){
      cat(filename, "has SNP in it! ")
      
      #save old exp file with name "AT1G01010_SNPin.txt"
      write.table(rawexp, file = paste0(location_snpCor, gsub(".txt", "_SNPin.txt", filename)))
      
      #save new exp file with name "AT1G01010.txt"
      rawexp <- rawexp[-which(rownames(rawexp) %in% snp), ]
      write.table(rawexp, file = paste0(location_snpCor, filename))
      
      cat(" old exp, new exp files saved! ")
      
      if(any(grepl(gsub(".txt", "", filename), dir(location_FM)))){
        #save new env file with name "AT1G01010_FM_Env.txt"
        env <- read.table(paste0(location_FM, gsub(".txt", "_FM_Env.txt", filename)), row.names=1, header=T)
        env <- env[-which(rownames(env) %in% snp), ]
        write.table(env, file = paste0(location_snpCor, gsub(".txt", "_FM_Env.txt", filename)))
        
        #save new qtl file with name "AT1G01010_FM_QTL.txt"
        qtl <- read.table(paste0(location_FM, gsub(".txt", "_FM_QTL.txt", filename)), row.names=1, header=T)
        qtl <- qtl[-which(rownames(qtl) %in% snp), ]
        write.table(qtl, file = paste0(location_snpCor, gsub(".txt", "_FM_QTL.txt", filename)))
        
        #save new int file with name "AT1G01010_FM_Int.txt"
        int <- read.table(paste0(location_FM, gsub(".txt", "_FM_Int.txt", filename)), row.names=1, header=T)
        int <- int[-which(rownames(int) %in% snp), ]
        write.table(int, file = paste0(location_snpCor, gsub(".txt", "_FM_Int.txt", filename)))
        
        cat(" new env, new qtl and new int files saved!")
      }
      cat("\n")
    } else{
      cat(filename, "no SNP\n")
    }
  }
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s!\n\n")
}


#*************************************************************** check part ***************************************************************#
filename="AT5G01400.txt"
rawexp <- read.table(paste0(location, filename), row.names=1, header=T)
dim(rawexp)
rawSNP <- read.table(paste0(location_snpCor, gsub(".txt", "_SNPin.txt", filename)), row.names=1, header=T)
dim(rawSNP)
newraw <- read.table(paste0(location_snpCor, filename), row.names=1, header=T)
dim(newraw)
newenv <- read.table(paste0(location_snpCor, gsub(".txt", "_FM_Env.txt", filename)), row.names=1, header=T)
dim(newenv)
any(rownames(rawexp) %in% snp)
any(rownames(rawSNP) %in% snp)
any(rownames(newraw) %in% snp)
any(rownames(newenv) %in% snp)


#************************************************************ move and replace ************************************************************#
for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  location_FM <- paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/")
  location_snpCor <- paste0("Data/snpCorrection/chr", chr, "_norm_hf_cor/")

  #first, move "AT1G01010_SNPin.txt"
  snps <- dir(location_snpCor)[grepl("_SNPin.txt", dir(location_snpCor))]
  file.rename(paste0(location_snpCor, snps), paste0(location, snps))

  #next, move "AT1G01010_FM_Env.txt", "AT1G01010_FM_QTL.txt" and "AT1G01010_FM_Int.txt"
  fms <- dir(location_snpCor)[grepl("FM", dir(location_snpCor))]
  file.rename(paste0(location_snpCor, fms), paste0(location_FM, fms))

  #finally, only "AT1G01010.txt" left, and move them
  file.rename(paste0(location_snpCor, dir(location_snpCor)), paste0(location, dir(location_snpCor)))
}

#********************************************************** delete no probe files **********************************************************#
for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  location_FM <- paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/")
  
  #filename = "AT1G01010.txt"
  for(filename in dir(location)[grepl(".txt", dir(location)) & !grepl("_QTL.txt", dir(location)) & !grepl("_SNPin.txt", dir(location))]){
    rawexp <- read.table(paste0(location, filename), row.names=1, header=T)
    if(nrow(rawexp) == 0){
      cat(filename, "no probe left!\n")
      
      #remove no probe gene files from exp folder
      file.rename(paste0(location, filename), paste0("no probe genes/", filename))
      
      #remove no probe gene files from mapping folder. If current gene have no mapping files, it won't move any files!
      file.rename(paste0(location_FM, dir(location_FM)[grepl(gsub(".txt", "", filename), dir(location_FM))]), paste0("no probe genes/", dir(location_FM)[grepl(gsub(".txt", "", filename), dir(location_FM))]))
    }
  }
  
}

#********************************************************** calculate gene numbers **********************************************************#
ngene <- NULL
nexp <- NULL
for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  location_FM <- paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/")
  
  ngene <- c(ngene, length(dir(location)[grepl(".txt", dir(location)) & !grepl("_QTL.txt", dir(location)) & !grepl("_SNPin.txt", dir(location))]))
  nexp <- c(nexp, length(dir(location_FM)[grepl("QTL", dir(location_FM))]))
}
