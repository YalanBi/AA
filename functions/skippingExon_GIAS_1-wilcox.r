#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 02-07-2013
# first written: 21-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#**************************************** this is the final version for testing GENETIC/INTERACTION regulated skipping exons! ^_^ ****************************************#
#****************************************************** skipping exon: cassette exons + the first/last spliced exon ******************************************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#

setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]
#load genotype file
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)
#load exp genes
load(file="Data/ExpGenes/expGenes_final.Rdata")

#direction selection
probesDir <- function(exp_data = rawexp){
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
#            whichTest <- wilcox.test/ t.test; one-side!!!
testDffBtwParts <- function(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind, whichTest=wilcox.test, alternative="less", verbose=FALSE){
  testPart <- as.numeric(unlist(exp_data[testProbes, ind]))
  if(verbose) cat("We are testProbes:", testProbes, "\n")
  restPart <- as.numeric(unlist(exp_data[restProbes, ind]))
  if(verbose) cat("We are restProbes:", restProbes, "\n")
  return(-log10(whichTest(testPart, restPart, alternative) $ p.value))
}
#testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind, whichTest=wilcox.test, alternative="less", verbose=TRUE)

#skipping exon test
#annotation: "goal" could be "skippingExon", "cassetteExon" and "skipping53Exon"
#            "toTest" could be "QTL" and "Int"
splicingTestSE_AS <- function(filename, goal="skippingExon", toTest="QTL", threTest=8, verbose=FALSE, ...){
  if(verbose) cat("now is testing", goal, "!\n")
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")
  
  #for 5'/3' AS, at least 2 exons in a gene!!!
  #NOTE: whether have length(exonID) >= 6 or not, it won't change the number of test, so that no change to Bonferroni crection and threshold
  #      having length(exonID) >= 6 is just for saving time, skip loading eQTL/Int for genes having 2 exons but in 1 exon there are only 1/2 probes
  if(length(uniqueExon) >= 2 && length(exonID) >= 6){
    if(verbose) cat(filename, "have", length(uniqueExon), "exons, and", length(exonID), "exonProbes!\n")
    
    testQTL <- read.table(paste0("Data/FullModel/chr", chr, "_norm_hf_cor_FM/", filename, "_FM_", toTest, ".txt"), row.names=1, header=TRUE)
    
    resmatrix <- NULL
    rownameList <- NULL
    
    if(goal == "skippingExon") testExonRange <- uniqueExon
    if(goal == "cassetteExon") testExonRange <- uniqueExon[-c(1, length(uniqueExon))]
    if(goal == "skipping53Exon") testExonRange <- uniqueExon[c(1, length(uniqueExon))]
    
    for(testExon in testExonRange){
      ind <- exonID[rawexp[exonID, "tu"] == testExon]
      #if(verbose) cat(testExon, "has probes", ind, "\n")
      
      testProbes <- ind
      #if(verbose) cat("I'm testpart, I have p", testProbes, "\n")
      restProbes <- exonID[!exonID %in% ind]
      #if(verbose) cat("I'm restpart, I have p", restProbes, "\n")
      
      #at least 3 probes in each group, try to avoid this case---one is highly expressed, the other is lowly expressed...
      if(length(testProbes) >= 3 && length(restProbes) >= 3 && any(apply(testQTL[testProbes, ] >= threTest, 2, sum) >= 3)){
        if(verbose) cat("\t***I'm", testExon, "sepProbe is p", max(testProbes), ", have >= 3 good probes in each group, and I have sig", toTest, ", ready for test!\n")
        
        m <- which(apply(testQTL[testProbes, ] >= threTest, 2, sum) >= 3)[which.max(apply(as.matrix(testQTL[testProbes, which(apply(testQTL[testProbes, ] >= threTest, 2, sum) >= 3)]), 2, sum))]
        if(verbose) cat(filename, testExon, "the most sig marker is", m, "among possible ones", which(apply(testQTL[testProbes, ] >= threTest, 2, sum) >= 3), "\n")
        geno1 <- which(geno[ ,m] == 1)
        geno2 <- which(geno[ ,m] == 2)
        
        #>= 3 probes left in each group and have sig eQTL on testExon, remember the gene name and do t.test
        rownameList <- c(rownameList, filename)
        res <- c(max(testProbes), as.numeric(m))
        for(env in 1:4){
          ind_env <- which(as.numeric(menvironment) == env)
          envGroup1 <- ind_env[ind_env %in% geno1]
          envGroup2 <- ind_env[ind_env %in% geno2]
          if(length(envGroup1) > 0 && length(envGroup2) > 0){
            res <- c(res, testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=envGroup1, verbose, ...), testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=envGroup2, verbose, ...))
            #if(verbose) cat("env", env, ": gt1 -", testDffBtwParts(exp_data = rawexp[ ,17:164], ind = envGroup1, ...), "; gt2 -", testDffBtwParts(exp_data = rawexp[ ,17:164], ind = envGroup2, ...),"\n")
          } else{
            if(verbose) cat("in env", env, ", one genotype have no RILs\n")
            res <- c(res, -1, -1)
          }
        }
        resmatrix <- rbind(resmatrix, res)
      } else if(verbose) cat("\t***I'm", testExon, ", not enough probes in each part or have no sig QTL on me T^T\n")
    }
    rownames(resmatrix) <- rownameList
    return(resmatrix)
  } else if(verbose) cat(filename, "we don't have enough exons or less than 6 exon probes T^T\n")
}
#splicingTestSE_AS(filename, goal="skippingExon", toTest="QTL", threTest=8, whichTest=wilcox.test, alternative="less", verbose=TRUE)


#test SE_GAS for chr 1-5
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- splicingTestSE_AS(filename, goal="skippingExon", toTest="QTL", threTest=8, whichTest=wilcox.test, alternative="less")
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      cat(filename, "is tested\n")
    }
  }
  if(!is.null(resmatrix)){
    colnames(resmatrix) <- c("sepProbe", "sigMarker", "6H/gt1", "6H/gt2", "Dry_AR/gt1", "Dry_AR/gt2", "Dry_Fresh/gt1", "Dry_Fresh/gt2", "RP/gt1", "RP/gt2")
    write.table(resmatrix, file=paste0("Data/geneticsAS/SE_chr", chr, "_GAS_wt_less.txt"), sep="\t") #********** change!!! **********#
  } else cat("\tNO TEST!\n")
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}

#test SE_IAS for chr 1-5
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- splicingTestSE_AS(filename, goal="skippingExon", toTest="Int", threTest=11.6, whichTest=wilcox.test, alternative="less")
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      cat(filename, "is tested\n")
    }
  }
  if(!is.null(resmatrix)){
    colnames(resmatrix) <- c("sepProbe", "sigMarker", "6H/gt1", "6H/gt2", "Dry_AR/gt1", "Dry_AR/gt2", "Dry_Fresh/gt1", "Dry_Fresh/gt2", "RP/gt1", "RP/gt2")
    write.table(resmatrix, file=paste0("Data/geneticsAS/SE_chr", chr, "_IAS_wt_less.txt"), sep="\t") #********** change!!! **********#
  } else cat("\tNO TEST!\n")
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}
