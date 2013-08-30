#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 23-07-2013
# first written: 23-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#******************************************************** this is the final version for testing skipping exons ^_^ *******************************************************#
#****************************************************** skipping exon: cassette exons + the first/last spliced exon ******************************************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#
#main idea: expGenes >= 2 exons
#           each exons >= 2 probes
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
testDffBtwParts <- function(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind, alternative, verbose=FALSE){
  testPart <- as.numeric(unlist(exp_data[testProbes, ind]))
  if(verbose) cat("We are testProbes:", testProbes, "\n")
  restPart <- as.numeric(unlist(exp_data[restProbes, ind]))
  if(verbose) cat("We are restProbes:", restProbes, "\n")
  return(-log10(wilcox.test(testPart, restPart, alternative) $ p.value))
}
#testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env, alternative, verbose=TRUE)

#skipping exon test
splicingTestSE <- function(filename, P=3, verbose=FALSE){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")
  intronID <- probes_dir[grepl("intron", rawexp[probes_dir,"tu"])]
  unExpInPrb <- intronID[apply(rawexp[intronID, 17:164], 1, median) < 5]
  expExPrb <- exonID[apply(rawexp[exonID, 17:164], 1, median) >= 5]
  
  #for skipping exons, at least 2 exons in a gene!!!
  if(length(uniqueExon) >= 2 && length(exonID) >= 2*P && length(expExPrb) >= P && length(unExpInPrb) >= P){
    if(verbose) cat(filename, "have", length(uniqueExon), "exons,", length(exonID), "exPrbs, Exp exPrbs(", expExPrb, ") and notExp inPrbs(", unExpInPrb, ")!\n")
    
    testRange <- NULL
    for(tu in uniqueExon){
      ind <- exonID[rawexp[exonID, "tu"] == tu]
      if(length(ind) > 0){
        if(median(unlist(rawexp[ind, 17:164])) < 5 && length(ind) >= P){
          if(verbose) cat("\tkeep", tu, "and put it into testPart: median =", median(unlist(rawexp[ind, 17:164])), "and having", length(ind), "probes\n")
          testRange <- c(testRange, tu)
        } else if(verbose) cat("\tbyebye", tu, ", can't be a testExon\n")
      } else if(verbose) cat("\tignore", tu, ", NO right_dir probes\n")
    }
    if(!is.null(testRange)){
      if(verbose) cat("we have exons(", testRange, ") to be tested, median of which is < 5, continue!\n")
      
      resmatrix <- NULL
      rownameList <- NULL
      
      for(toTest in testRange){
        testProbes <- exonID[rawexp[exonID, "tu"] == toTest]
        restExpExP <- expExPrb[!expExPrb %in% testProbes]
        
        if(length(restExpExP) >= P){
          if(verbose) cat("\t***I'm", toTest, ", sepProbe is p", max(testProbes), ", enough testProbes(", testProbes, ") and restExpExonProbes(", restExpExP, "), test me!\n")
          #>= P probes left in each group, remember the gene name and do wilcox.test
          rownameList <- c(rownameList, filename)
          res <- max(testProbes)
          for(env in 1:4){
            ind_env <- which(as.numeric(menvironment) == env)
            res <- c(res, testDffBtwParts(rawexp[ ,17:164], testProbes, unExpInPrb, ind_env, "greater"), testDffBtwParts(rawexp[ ,17:164], testProbes, restExpExP, ind_env, "less"))
            #if(verbose) cat("env", env, ":", testDffBtwParts(rawexp[ ,17:164], testProbes, unExpInPrb, ind_env, "greater"), testDffBtwParts(rawexp[ ,17:164], testProbes, restExpExP, ind_env, "less"), "\n")
          }
          resmatrix <- rbind(resmatrix, res)
        } else if(verbose) cat("\t***I'm", toTest, "not enough probes in restPart, skip me T^T\n")
      }
      rownames(resmatrix) <- rownameList
      return(resmatrix)
    } else if(verbose) cat("\twe don't have any exons to be tested T^T\n")
  } else if(verbose) cat(filename, "don't have enough exons/exPrbs/Exp exPrbs/notExp inPrbs T^T\n")
}
#splicingTestSE(filename, P=2, verbose=TRUE)


#test for chr 1-5
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- splicingTestSE(filename, P=3)
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      cat(filename, "is tested\n")
    }
  }
  if(!is.null(resmatrix)){
    colnames(resmatrix) <- c("sepProbe", "6H/In", "6H/Ex", "Dry_AR/In", "Dry_AR/Ex", "Dry_Fresh/In", "Dry_Fresh/Ex", "RP/In", "RP/Ex")
    write.table(resmatrix, file=paste0("Data/AS/SE_chr", chr, "_Ex+In.txt"), sep="\t") #***************************** change with goal !!! *****************************#
  } else cat("\tNO TEST!\n")
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}
