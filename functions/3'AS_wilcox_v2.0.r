#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 29-07-2013
# first written: 22-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#********************************************************* this is the final version for testing AS at 3 site ^_^ ********************************************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#
#main idea:
#minimum of 2*P probes in this exon; minimum of (2*P+1) probes in this gene
#Test how to split (using highest difference between two groups)
#Split into two groups
#test if every group has P probes
#YES -> T-Test the test group (5" start) against (the rest part of first exon + all the probes from the other expExons in the gene)
#NO -> continue

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

#to find max difference between each probe in exp, for grouping
findSepPoint <- function(toGroup=ind, exp_data=rawexp[ ,17:164], verbose=FALSE){
  dffs <- NULL
  for(n in 2:length(toGroup)){
    dff <- sum(apply(exp_data[toGroup[length(toGroup):n], ], 2, median)-apply(exp_data[toGroup[1:(n-1)], ], 2, median))
    dffs <- c(dffs, dff)
    if(verbose) cat("difference between probe(", toGroup[length(toGroup):n], ") and probe(", toGroup[1:(n-1)], ")is", dff, "\n")
  }
  if(min(dffs) < 0){
    if(verbose) cat("so dffList is:", dffs, "\n", "and edge probes are p", toGroup[which.min(dffs)], "and p", toGroup[which.min(dffs)+1], ", dff =", min(dffs), "\n")
    return(which.min(dffs))
  } else return() #we want the testPart is lower than the other part of this exon, otherwise it is decay/decrease
}
#findSepPoint(toGroup=ind, exp_data=rawexp[ ,17:164], verbose=TRUE)

#test the difference between 2 groups
#annotation: unlist -> use all individuals to do test, better than mean/median
#            whichTest <- wilcox.test/ t.test; one-side!!!
testDffBtwParts <- function(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind, verbose=FALSE){
  testPart <- as.numeric(unlist(exp_data[testProbes, ind]))
  if(verbose) cat("We are testProbes:", testProbes, "\n")
  restPart <- as.numeric(unlist(exp_data[restProbes, ind]))
  if(verbose) cat("We are restProbes:", restProbes, "\n")
  return(-log10(wilcox.test(testPart, restPart, alternative="less") $ p.value))
}
#testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env, verbose=TRUE)

#3' AS test
splicingTestAt3 <- function(filename, P=2, verbose=FALSE){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")
  
  restProbes <- NULL
  for(tu in uniqueExon[-length(uniqueExon)]){
    ind <- exonID[rawexp[exonID, "tu"] == tu]
    if(length(ind) > 0 && median(unlist(rawexp[ind, 17:164])) >= 5){
      if(verbose) cat("\tput", tu, "into restPart: median =", median(unlist(rawexp[ind, 17:164])), ">= 5\n")
      restProbes <- c(restProbes, ind)
    } else if(verbose) cat("\tignore", tu, ", NO right_dir probes/median < 5\n")
  }
  
  ind <- sort(exonID[rawexp[exonID, "tu"] == uniqueExon[length(uniqueExon)]], decreasing=TRUE)
  #if(verbose) cat(uniqueExon[length(uniqueExon)], "has probes", ind, "\n")
  
  #for 5'/3' AS, at least 2 exons in a gene!!! And at least 2*P probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
  if(length(uniqueExon) >= 2 && length(exonID) >= (2*P+1) && length(ind) >= 2*P && !is.null(restProbes)){
    if(verbose) cat(filename, "have", length(uniqueExon), "exons,", length(exonID), "exonProbes,", length(ind), "exonProbes in lastExon and have expExons, can be tested for 3'AS!\n")
    
    sepPoint <- findSepPoint(toGroup=ind, exp_data=rawexp[ ,17:164], verbose=FALSE)
    if(!is.null(sepPoint) && length(ind[1:sepPoint]) >= P && length(ind[-(1:sepPoint)]) >= P && median(unlist(rawexp[ind[sepPoint:1], 17:164])) >= 5 && median(unlist(rawexp[ind[-(1:sepPoint)], 17:164])) < 5){
      if(verbose) cat(" =>separate between p", ind[sepPoint], "and p", ind[sepPoint+1], ":",
      "partI(", ind[-(1:sepPoint)], ") median is", median(unlist(rawexp[ind[-(1:sepPoint)], 17:164])), "; partII(", ind[sepPoint:1], ") median is", median(unlist(rawexp[ind[sepPoint:1], 17:164])),
      "\nready for test!\n")
      
      testProbes <- ind[-(1:sepPoint)]
      if(verbose) cat("I'm testpart, I have p", testProbes, "\n")
      restProbes <- c(restProbes, ind[sepPoint:1])
      if(verbose) cat("I'm restpart, I have p", restProbes, "\n")
      
      #NOTE: min(ind[sepPoint], ind[sepPoint+1]) <- the probe just before the gap, for making plot.
      #      in 5'AS, it is the last probe of higher part in the first exon; in 3'AS, it is the last probe of the lower part in the last exon
      res <- ind[sepPoint+1]
      #>= P probes left in each group of the first exon, remember the gene name and do wilcox.test
      for(env in 1:4){
        ind_env <- which(as.numeric(menvironment) == env)
        res <- c(res, testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env))
        #if(verbose) cat("env", env, ":", testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env), "\n")
      }
      return(res)
    } else if(verbose) cat(" =>last exon decreases/not enough exonProbes in each part after separation T^T\n")
  } else if(verbose) cat("we don't have enough exons/exonProbes/exonProbes in lastExon/no expExon T^T\n")
}
#splicingTestAt3(filename, P=2, verbose=TRUE)


#test 3'AS for chr 1-5
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  rownameList <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- splicingTestAt3(filename, P=2)
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      rownameList <- c(rownameList, filename)
      cat(filename, "is tested\n")
    }
  }
  if(!is.null(resmatrix)){
    rownames(resmatrix) <- rownameList
    colnames(resmatrix) <- c("sepProbe", "6H", "Dry_AR", "Dry_Fresh", "RP")
    write.table(resmatrix, file=paste0("Data/AS/3'AS_chr", chr, "_wt_less_p2.txt"), sep="\t") #**************************** change with goal!!! ****************************#
  } else cat("\tNO TEST!\n")
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}
