#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 22-07-2013
# first written: 22-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#******************************************************** this is the final version for testing consistent Int ^_^ *******************************************************#

setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]
#load map file---markers' Morgan position
mpos <- read.table("refined map/map.physical.c15o15f.txt", row.names=1, header=T)
#load exp genes
load(file="Data/ExpGenes/expGenes_final.Rdata")

#test that each tu is cis-/trans- eQTL, or nothing...
consistent_Int <- function(chr, fnInt, intThre=11.6, fnPos, verbose=FALSE){
  res <- NULL
  for(e in 1:nrow(fnInt)){
    topM <- which.max(fnInt[e, ])
    if(verbose) cat("I'm m", topM, ", top marker of", rownames(fnInt)[e], ", Int =", fnInt[e, topM], "!\n")
    if(fnInt[e, topM] >= intThre){
      if(verbose) cat("\tsig Int and continue!\n")
      
      if(mpos[topM, "chr"] != chr){
        res <- c(res, "trans")
        if(verbose) cat("\ttop marker and gene are on dff chr, so TRANS!\n")
      } else if(abs(mpos[topM, "bp"]-fnPos[e]) > 5000000){
        res <- c(res, "trans")
        if(verbose) cat("\ttop marker and gene are on same chr but dstnc is > 5Mb, so TRANS!\n")
      } else{
        res <- c(res, "cis")
        if(verbose) cat("\ttop marker and gene are on same chr and dstnc is <= 5Mb, so CIS!\n")
      }
    } else{
      if(verbose) cat("\tnot sig Int T^T, quit...\n")
      return(NULL)
    }
  }
  if(all(res == "cis")) return("cis")
  else if(all(res == "trans")) return("trans")
  else return(NULL)
}


#get consistent Int genes
consistentInt <- list()
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  rawexp <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_1tu1probe.txt"), row.names=1, header=TRUE)
  int <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_Int.txt"), row.names=1, header=TRUE)
  
  genenames <- expGeneList[[chr]]
  for(filename in genenames){
    fnPos <- rawexp[grep(filename, rownames(rawexp)), 1]
    fnInt <- int[grep(filename, rownames(int)), ]
    
    res <- consistent_Int(chr, fnInt, intThre=11.6, fnPos)
    if(!is.null(res)){
      if(res == "cis") consistentInt$cis <- c(consistentInt$cis, filename)
      if(res == "trans") consistentInt$trans <- c(consistentInt$trans, filename)
      cat(filename, "found have consistent", res, "eQTL!\n")
    }
  }
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s!\n")
}
save(consistentInt, file=paste0("Data/testQTL/consistent_Int.Rdata"))
