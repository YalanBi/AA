#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 22-07-2013
# first written: 17-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#*********************************************************** this is the final version for testing main Int ^_^ **********************************************************#

setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]
#load map file---markers' Morgan position
mpos <- read.table("refined map/map.physical.c15o15f.txt", row.names=1, header=T)
#load exp genes
load(file="Data/ExpGenes/expGenes_final.Rdata")

#test that each tu is cis-/trans- Int, or nothing...
main_Int <- function(chr, fnInt, intThre=11.6, fnPos, verbose=FALSE){
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
      res <- c(res, "")
      if(verbose) cat("\tnot sig Int T^T, next one...\n")
    }
  }
  if("cis" %in% res && "trans" %in% res) return("both")
  else if("cis" %in% res) return("cis")
  else if("trans" %in% res) return("trans")
  else return(NULL)
}


#get main Int genes
mainInt <- list()
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  rawexp <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_1tu1probe.txt"), row.names=1, header=TRUE)
  int <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_Int.txt"), row.names=1, header=TRUE)
  
  genenames <- expGeneList[[chr]]
  for(filename in genenames){
    fnPos <- rawexp[grep(filename, rownames(rawexp)), 1]
    fnInt <- int[grep(filename, rownames(int)), ]
    
    res <- main_Int(chr, fnInt, intThre=11.6, fnPos)
    if(!is.null(res)){
      if(res == "cis") mainInt$cis <- c(mainInt$cis, filename)
      if(res == "trans") mainInt$trans <- c(mainInt$trans, filename)
      if(res == "both") mainInt$both <- c(mainInt$both, filename)
      cat(filename, "found have", res, "Int!\n")
    }
  }
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s!\n")
}
save(mainInt, file=paste0("Data/testQTL/main_Int.Rdata"))
