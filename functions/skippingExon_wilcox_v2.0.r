#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 18-07-2013
# first written: 18-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#******************************************************** this is the final version for testing skipping exons ^_^ *******************************************************#
#****************************************************** skipping exon: cassette exons + the first/last spliced exon ******************************************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#
#main idea: expGenes >= 2 exons
#           each exons >= 3 probes
#           testExon: meadian < 5; restExons: median >= 5; wilcox.test

setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]

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

#test the difference between 2 groups
#annotation: unlist -> use all individuals to do test, better than mean/median
testDffBtwParts <- function(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind, alternative="less", verbose=FALSE){
  testPart <- as.numeric(unlist(exp_data[testProbes, ind]))
  if(verbose) cat("We are testProbes:", testProbes, "\n")
  restPart <- as.numeric(unlist(exp_data[restProbes, ind]))
  if(verbose) cat("We are restProbes:", restProbes, "\n")
  return(-log10(wilcox.test(testPart, restPart, alternative) $ p.value))
}
#testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env, alternative="less", verbose=TRUE)

#skipping exon test
splicingTestSE <- function(filename, P=2, verbose=FALSE, ...){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")
  
  #for skipping exons, at least 2 exons in a gene!!!
  if(length(uniqueExon) >= 2 && length(exonID) >= 2*P){
    if(verbose) cat(filename, "have", length(uniqueExon), "exons, and", length(exonID), "exonProbes!\n")
    
    restProbes <- NULL
    testRange <- NULL
    for(tu in uniqueExon){
      ind <- exonID[rawexp[exonID, "tu"] == tu]
      if(length(ind) > 0){
        if(median(unlist(rawexp[exonID[rawexp[exonID, "tu"] == tu], 17:164])) >= 5){
          if(verbose) cat("\tkeep", tu, "and put it into restPart\n")
          restProbes <- c(restProbes, ind)
        } else if(length(ind) >= P){
          if(verbose) cat("\tkeep", tu, "and put it into testPart\n")
          testRange <- c(testRange, tu)
        } else if(verbose) cat("\tbyebye", tu, "forever, too few probes to be a testExon\n")
      } else if(verbose) cat("\tbyebye", tu, "forever, NO right_dir probes\n")
    }
    if(!is.null(testRange) && length(restProbes) >= P){
      if(verbose) cat("we have exons(", testRange, ") to be tested, median of which is < 5, and enough restProbes", restProbes, "!\n")
      
      resmatrix <- NULL
      rownameList <- NULL
      
      for(toTest in testRange){
        testProbes <- exonID[rawexp[exonID, "tu"] == toTest]
        if(verbose) cat("\t***I'm", toTest, "sepProbe is p", max(testProbes), ", have >=", P, "good probes in each group, ready for test!\n")
        #>= P probes left in each group, remember the gene name and do wilcox.test
        rownameList <- c(rownameList, filename)
        res <- max(testProbes)
        for(env in 1:4){
          ind_env <- which(as.numeric(menvironment) == env)
          res <- c(res, testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env, verbose, ...))
          #if(verbose) cat("env", env, ":", testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env, ...), "\n")
        }
        resmatrix <- rbind(resmatrix, res)
      }
      rownames(resmatrix) <- rownameList
      return(resmatrix)
    } else if(verbose) cat("\twe don't have exons to be tested, or not enough restProbes T^T\n")
  } else if(verbose) cat("we don't have enough exons T^T\n")
}
#splicingTestSE(filename, P=2, alternative="less", verbose=TRUE)


#test for chr 1-5
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- splicingTestSE(filename, P=2, alternative="less")
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      cat(filename, "is tested\n")
    }
  }
  if(!is.null(resmatrix)){
    colnames(resmatrix) <- c("sepProbe", "6H", "Dry_AR", "Dry_Fresh", "RP")
    write.table(resmatrix, file=paste0("Data/AS/SE_chr", chr, "_wt_less.txt"), sep="\t") #***************************** change with goal !!! *****************************#
  } else cat("\tNO TEST!\n")
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}
