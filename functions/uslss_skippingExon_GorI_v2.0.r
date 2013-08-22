#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 25-07-2013
# first written: 22-07-2013
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
#load genotype file
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)
#load testable genes
load(file="Data/AS/SE_testableGenes.Rdata")

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
testDffBtwParts <- function(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind, verbose=FALSE){
  testPart <- as.numeric(unlist(exp_data[testProbes, ind]))
  if(verbose) cat("We are testProbes:", testProbes, "\n")
  restPart <- as.numeric(unlist(exp_data[restProbes, ind]))
  if(verbose) cat("We are restProbes:", restProbes, "\n")
  return(-log10(wilcox.test(testPart, restPart, alternative="less") $ p.value))
}
#testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env, verbose=TRUE)

#G/I regulated skipping exon test
#annotation: "goal" could be "QTL" and "Int"
splicingGorISE <- function(filename, goal="QTL", threTest=8, fnPart, P=3, verbose=FALSE){
  if(verbose) cat("now is testing", filename, "for", goal, "regulated skipping exons!\n")
  
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")
  
  #for skipping exons, at least 2 exons in a gene!!!
  #NOTE: whether have length(exonID) >= 2*P or not, it won't change the number of test, so that no change to Bonferroni crection and threshold
  #      having length(exonID) >= 2*P is just for saving time, skip loading eQTL/Int for genes having 2 exons but in 1 exon there are only 1/2 probes
#  if(length(uniqueExon) >= 2 && length(exonID) >= 2*P){
#    if(verbose) cat(filename, "have", length(uniqueExon), "exons, and", length(exonID), "exonProbes!\n")
    
    testRange <- NULL
    restProbes <- NULL
    for(tu in uniqueExon){
      ind <- exonID[rawexp[exonID, "tu"] == tu]
      if(length(ind) > 0){
        if(median(unlist(rawexp[ind, 17:164])) >= 5){
          if(verbose) cat("\tput", tu, "into restPart: median =", median(unlist(rawexp[ind, 17:164])), "> 5\n")
          restProbes <- c(restProbes, ind)
        } else if(length(ind) >= P){
          if(verbose) cat("\tkeep", tu, "and put it into testPart: median =", median(unlist(rawexp[ind, 17:164])), "and having", length(ind), "probes\n")
          testRange <- c(testRange, tu)
        } else if(verbose) cat("\tbyebye", tu, ", not enough probes to be a testExon\n")
      } else if(verbose) cat("\tignore", tu, ", NO right_dir probes\n")
    }
    if(!is.null(testRange) && length(restProbes) >= P){
      if(verbose) cat("we have exons(", testRange, ") to be tested and enough restProbes(", restProbes, "), continue testing having sig", goal, "or not!\n")
      
      resmatrix <- NULL
      rownameList <- NULL
      
      for(toTest in testRange){
        testProbes <- exonID[rawexp[exonID, "tu"] == toTest]
        
        if(max(fnPart[paste0(filename, "_", toTest), ]) >= threTest){
          if(verbose) cat("\t***I'm", toTest, ", enough testProbes(", testProbes, ") and restProbes(", restProbes, "), and sig", goal, ", test me!\n")
          
          m <- which.max(fnPart[paste0(filename, "_", toTest), ])
          if(verbose) cat(filename, toTest, "the most sig marker is", m, "\n")
          geno1 <- which(geno[ ,m] == 1)
          geno2 <- which(geno[ ,m] == 2)
          
          #>= P probes left in each group and have sig eQTL on testExon, remember the gene name and do wilcox.test
          rownameList <- c(rownameList, filename)
          res <- c(max(testProbes), m)
          for(env in 1:4){
            ind_env <- which(as.numeric(menvironment) == env)
            envGroup1 <- ind_env[ind_env %in% geno1]
            envGroup2 <- ind_env[ind_env %in% geno2]
            if(length(envGroup1) > 0 && length(envGroup2) > 0){
              res <- c(res, testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=envGroup1), testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=envGroup2))
              #if(verbose) cat("env", env, ": gt1 -", testDffBtwParts(exp_data=rawexp[ ,17:164], ind=envGroup1), "; gt2 -", testDffBtwParts(exp_data=rawexp[ ,17:164], ind=envGroup2), "\n")
            } else{
              if(verbose) cat("in env", env, ", one genotype have no RILs\n")
              res <- c(res, -1, -1)
            }
          }
          resmatrix <- rbind(resmatrix, res)
        } else if(verbose) cat("\t***I'm", toTest, "no sig", goal, "on me, skip me T^T\n")
      }
      rownames(resmatrix) <- rownameList
      return(resmatrix)
    } else if(verbose) cat("\twe don't have any exons to be tested/not enough restProbes T^T\n")
#  } else if(verbose) cat(filename, "don't have enough exons/exonProbes T^T\n")
}
#splicingGorISE(filename, goal="QTL", threTest=8, fnPart, P=3, verbose=TRUE)


#test GSE for chr 1-5
goal="QTL"# "Int"
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts testing", goal, "regulated skipping exons...\n")
  
  qtl <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_", goal, ".txt"), row.names=1, header=TRUE)
  
  genenames <- grep(paste0("AT", chr), testableGenes, value=TRUE)
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    fnPart <- abs(qtl[grep(filename, rownames(qtl)), ])
    
    res <- splicingGorISE(filename, goal="QTL", threTest=8, fnPart, P=3)
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      cat(filename, "is tested\n")
    }
  }
  if(!is.null(resmatrix)){
    colnames(resmatrix) <- c("sepProbe", "sigMarker", "6H/gt1", "6H/gt2", "Dry_AR/gt1", "Dry_AR/gt2", "Dry_Fresh/gt1", "Dry_Fresh/gt2", "RP/gt1", "RP/gt2")
    write.table(resmatrix, file=paste0("Data/geneticsAS/GSE_chr", chr, "_wt_p3.txt"), sep="\t") #***************************** change with goal !!! *****************************#
  } else cat("\tNO TEST!\n")
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}

#test ISE for chr 1-5
goal="Int"
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts testing", goal, "regulated skipping exons...\n")
  
  int <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_", goal, ".txt"), row.names=1, header=TRUE)
    
  genenames <- grep(paste0("AT", chr), testableGenes, value=TRUE)
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    fnPart <- abs(int[grep(filename, rownames(int)), ])
    
    res <- splicingGorISE(filename, goal="Int", threTest=11.6, fnPart, P=3)
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      cat(filename, "is tested\n")
    }
  }
  if(!is.null(resmatrix)){
    colnames(resmatrix) <- c("sepProbe", "sigMarker", "6H/gt1", "6H/gt2", "Dry_AR/gt1", "Dry_AR/gt2", "Dry_Fresh/gt1", "Dry_Fresh/gt2", "RP/gt1", "RP/gt2")
    write.table(resmatrix, file=paste0("Data/geneticsAS/ISE_chr", chr, "_wt_p3.txt"), sep="\t") #***************************** change with goal !!! *****************************#
  } else cat("\tNO TEST!\n")
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}
