#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 17-07-2013
# first written: 17-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#***************************************************** this is the final version for testing main/consistent Int! ^_^ ****************************************************#

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

test_Int <- function(goal, geneInt, intThre=11.6, verbose=FALSE){
  resInt <- colSums(abs(geneInt) >= intThre)
  
  if(goal == "main") P=1
  if(goal == "consistent") P=nrow(geneQTL)
  
  for(m in 1:716){
    if(resInt[m] >= P){
      if(verbose) cat("I'm", filename, "have", goal, "Int at marker", m, "^_^\n")
      return(TRUE)
    }
  }
  return(FALSE)
}

#get main Int genes
mainInt <- list()
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  int <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_Int.txt"), row.names=1, header=TRUE)
  
  genenames <- expGeneList[[chr]]
  for(filename in genenames){
    geneInt <- int[grep(filename, rownames(int)), ]
    if(test_Int(goal="main", geneInt, intThre=11.6)){
      mainInt[[paste0("chr", chr)]] <- c(mainInt[[paste0("chr", chr)]], filename)
      cat(filename, "found have main Int!\n")
    }
  }
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s!\n")
}
save(mainInt, file=paste0("Data/testQTL/main_Int.Rdata"))

#get consistent Int genes
consistentInt <- list()
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  int <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_Int.txt"), row.names=1, header=TRUE)
  
  genenames <- expGeneList[[chr]]
  for(filename in genenames){
    geneInt <- int[grep(filename, rownames(int)), ]
    if(test_Int(goal="consistent", geneInt, intThre=11.6)){
      consistentInt[[paste0("chr", chr)]] <- c(consistentInt[[paste0("chr", chr)]], filename)
      cat(filename, "found have consistent Int!\n")
    }
  }
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s!\n")
}
save(consistentInt, file=paste0("Data/testQTL/consistent_Int.Rdata"))
