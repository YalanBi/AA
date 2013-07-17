#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 17-07-2013
# first written: 17-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#***************************************************** this is the final version for testing main/consistent eQTL ^_^ ****************************************************#

setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]
#load exp genes
load(file="Data/ExpGenes/expGenes_final.Rdata")

#direction selection
probesDir <- function(exp_data=rawexp){
  if(unique(exp_data[ ,"strand"]) == "sense"){
    direction_id <- which(exp_data[ ,"direction"] == "reverse")
  }
  if(unique(exp_data[ ,"strand"]) == "complement"){
    direction_id <- which(exp_data[ ,"direction"] == "forward")
  }
  return(direction_id)
}

test_eQTL <- function(goal, geneQTL, qtlThre=8, geneInt, intThre=11.6, verbose=FALSE){
  resQTL <- colSums(abs(geneQTL) >= qtlThre)
  resInt <- colSums(abs(geneInt) >= intThre)
  
  if(goal == "main") P=1
  if(goal == "consistent") P=nrow(geneQTL)
  
  for(m in 1:716){
    if(resQTL[m] >= P && resInt[m] == 0){
      if(verbose) cat("I'm", filename, "have", goal, "eQTL at marker", m, "^_^\n")
      return(TRUE)
    }
  }
  return(FALSE)
}

#get main eQTL genes
mainQTL <- list()
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  qtl <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_QTL.txt"), row.names=1, header=TRUE)
  int <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_Int.txt"), row.names=1, header=TRUE)
  
  genenames <- expGeneList[[chr]]
  for(filename in genenames){
    geneQTL <- qtl[grep(filename, rownames(qtl)), ]
    geneInt <- int[grep(filename, rownames(int)), ]
    if(test_eQTL(goal="main", geneQTL, qtlThre=8, geneInt, intThre=11.6)){
      mainQTL[[paste0("chr", chr)]] <- c(mainQTL[[paste0("chr", chr)]], filename)
      cat(filename, "found have main eQTL!\n")
    }
  }
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s!\n")
}
save(mainQTL, file=paste0("Data/testQTL/main_eQTL.Rdata"))

#get consistent eQTL genes
consistentQTL <- list()
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  qtl <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_QTL.txt"), row.names=1, header=TRUE)
  int <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_Int.txt"), row.names=1, header=TRUE)
  
  genenames <- expGeneList[[chr]]
  for(filename in genenames){
    geneQTL <- qtl[grep(filename, rownames(qtl)), ]
    geneInt <- int[grep(filename, rownames(int)), ]
    if(test_eQTL(goal="consistent", geneQTL, qtlThre=8, geneInt, intThre=11.6)){
      consistentQTL[[paste0("chr", chr)]] <- c(consistentQTL[[paste0("chr", chr)]], filename)
      cat(filename, "found have consistent eQTL!\n")
    }
  }
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s!\n")
}
save(consistentQTL, file=paste0("Data/testQTL/consistent_eQTL.Rdata"))
