#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 10-07-2013
# first written: 08-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#******************************************************** this is the final version for testing skipping exons ^_^ *******************************************************#
#****************************************************** skipping exon: cassette exons + the first/last spliced exon ******************************************************#
#*********************************************************************** testing algorithm: ANOVA! ***********************************************************************#

setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]
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
#annotation: use ANOVA
testDffBtwParts <- function(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind, verbose=FALSE){
  testPart <- as.numeric(unlist(exp_data[testProbes, ind]))
  if(verbose) cat("We are testProbes:", length(testPart), "\n")
  restPart <- as.numeric(unlist(exp_data[restProbes, ind]))
  if(verbose) cat("We are restProbes:", length(restPart), "\n")
  
  response <- c(testPart, restPart)
  predictor <- c(rep(0,length(testPart)), rep(1, length(restPart)))
  model <- lm(response ~ as.numeric(predictor))
  if(model[[1]][2] < 0) return(0);
  -log10(anova(model)[[5]][1])
  #list(response, predictor)
}
#testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind, verbose=TRUE)

#skipping exon test
#annotation: "goal" could be "skippingExon", "cassetteExon" and "skipping53Exon"
splicingTestSE <- function(filename, goal="skippingExon", P=2, verbose=FALSE){
  if(verbose) cat("now is testing", goal, "!\n")
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
      
      #at least P probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
      if(length(testProbes) >= P && length(restProbes) >= P){
        if(verbose) cat("\t***I'm", testExon, "sepProbe is p", max(ind), ", have >=", P, "good probes in each group, ready for test!\n")
        
        #>= P probes left in each group, remember the gene name and do t.test
        rownameList <- c(rownameList, filename)
        res <- max(ind)
        for(env in 1:4){
          ind_env <- which(as.numeric(menvironment) == env)
          res <- c(res, testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env, verbose))
          #if(verbose) cat("env", env, ":", testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env), "\n")
        }
        resmatrix <- rbind(resmatrix, res)
      } else if(verbose) cat("\t***I'm", testExon, ", not enough probes in each part for test T^T\n")
    }
    rownames(resmatrix) <- rownameList
    return(resmatrix)
  } else if(verbose) cat("we don't have enough exons T^T\n")
}
#splicingTestSE(filename, goal="skippingExon", P=2, verbose=TRUE)


#test for chr 1-5
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- splicingTestSE(filename, goal="skippingExon", P=2)
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      cat(filename, "is tested\n")
    }
  }
  if(!is.null(resmatrix)){
    colnames(resmatrix) <- c("sepProbe", "6H", "Dry_AR", "Dry_Fresh", "RP")
    write.table(resmatrix, file=paste0("Data/AS/SE_chr", chr, "_ANOVA_less.txt"), sep="\t") #**************************** change with goal!!! ****************************#
  } else cat("\tNO TEST!\n")
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}
