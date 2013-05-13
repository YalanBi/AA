#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 02-05-2013
# first written: 02-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays/")



#*************************************************************** for genome ***************************************************************#
#count ngenes we have on each chr, after removing SNPs
nGenes <- NULL
for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  nGenes <- c(nGenes, length(dir(location)[grepl(".txt", dir(location)) & !grepl("SNP", dir(location)) & !grepl("QTL", dir(location))]))
}
nGenes
[1] 7683 4966 6249 4734 6896


#count ntus we have on each chr, after removing SNPs, for probes of right direction
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

nTus <- c(0, 0, 0, 0, 0)
for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  for(filename in dir(location)[grepl(".txt", dir(location)) & !grepl("SNP", dir(location)) & !grepl("QTL", dir(location))]){
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename), row.names=1, header=T)
    
    #get probes' and exons' ID
    probes_dir <- probesDir(rawexp)
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
    
    nTus[chr] <- nTus[chr] + length(unique(rawexp[exonID, "tu"]))
  }
}
nTus
[1] 34762 20142 25854 20294 30662


#count ntus we have on each chr, after removing SNPs, for probes of both directions
nTus <- c(0, 0, 0, 0, 0)
for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  for(filename in dir(location)[grepl(".txt", dir(location)) & !grepl("SNP", dir(location)) & !grepl("QTL", dir(location))]){
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename), row.names=1, header=T)
    
    nTus[chr] <- nTus[chr] + length(unique(rawexp[ ,"tu"][grepl("tu", rawexp[ ,"tu"])]))
  }
}
nTus
[1] 38181 22083 28318 22209 33489



#************************************************************** for ExpGenes **************************************************************#
#count nExpGenes we have on each chr, after removing SNPs
load(file="Data/fullModeMapping/expGenes.Rdata")
nExpGenes <- NULL
for(chr in 1:5){
  nExpGenes <- c(nExpGenes, length(expGeneList[[chr]]))
}
nExpGenes


#count nExpTus we have on each chr, after removing SNPs, for probes of right direction
load(file="Data/fullModeMapping/expGenes.Rdata")
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

nExpTus <- c(0, 0, 0, 0, 0)
for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  for(filename in expGeneList[[chr]]){
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename), row.names=1, header=T)
    
    #get probes' and exons' ID
    probes_dir <- probesDir(rawexp)
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
    
    nExpTus[chr] <- nExpTus[chr] + length(unique(rawexp[exonID, "tu"]))
  }
}
nExpTus
[1] 19087 11082 14830 11053 17443


#count nExpTus we have on each chr, after removing SNPs, for probes of both directions
load(file="Data/fullModeMapping/expGenes.Rdata")
nExpTus <- c(0, 0, 0, 0, 0)
for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  for(filename in expGeneList[[chr]]){
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename), row.names=1, header=T)
    
    nExpTus[chr] <- nExpTus[chr] + length(unique(rawexp[ ,"tu"][grepl("tu", rawexp[ ,"tu"])]))
  }
}
nExpTus
[1] 20886 12134 16228 12078 18984
