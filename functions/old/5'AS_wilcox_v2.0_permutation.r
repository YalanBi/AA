#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 31-07-2013
# first written: 31-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#********************************************************* this is the final version for testing AS at 5 site ^_^ ********************************************************#
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

#5' AS test
permutation5AS <- function(filename, P=2, verbose=FALSE){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")
  
  ind <- exonID[rawexp[exonID, "tu"] == uniqueExon[1]]
  #if(verbose) cat(uniqueExon[1], "has probes", ind, "\n")
  
  #for 5'/3' AS, at least 2 exons in a gene!!! And at least 2*P probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
  if(length(uniqueExon) >= 2 && length(exonID) >= (2*P+1) && length(ind) >= 2*P){
    if(verbose) cat(filename, "have", length(uniqueExon), "exons,", length(exonID), "exonProbes, and", length(ind), "exonProbes in firstExon, can be tested for 5'AS!\n")
    
    randomOrder <- sample(exonID, length(exonID), replace=FALSE)
    rawexp[exonID, "tu"] <- rawexp[randomOrder, "tu"]
    if(verbose) cat("new tu order:", as.character(rawexp[exonID, "tu"]), "\n")
    
    restProbes <- NULL
    for(tu in uniqueExon[-1]){
      tuPrbs <- exonID[rawexp[exonID, "tu"] == tu]
      if(length(tuPrbs) > 0 && median(unlist(rawexp[tuPrbs, 17:164])) >= 5){
        if(verbose) cat("\tput", tu, "into restPart: median =", median(unlist(rawexp[tuPrbs, 17:164])), ">= 5\n")
        restProbes <- c(restProbes, tuPrbs)
      } else if(verbose) cat("\tignore", tu, ", NO right_dir probes/median < 5\n")
    }
    
    if(!is.null(restProbes)){
      if(verbose) cat(" after permutation, have expExons and restProbes are", restProbes, "\n")
      
      ind <- exonID[rawexp[exonID, "tu"] == uniqueExon[1]] #new exonProbes in first exon
      sepPoint <- findSepPoint(toGroup=ind, exp_data=rawexp[ ,17:164], verbose=FALSE)
      if(!is.null(sepPoint) && sepPoint > 1 && sepPoint < (length(ind)-1) && median(unlist(rawexp[ind[1:sepPoint], 17:164])) >= 5 && median(unlist(rawexp[ind[-(1:sepPoint)], 17:164])) < 5){
        if(verbose) cat(" =>separate between p", ind[sepPoint], "and p", ind[sepPoint+1], ":",
        "partI(", ind[1:sepPoint], ") median is", median(unlist(rawexp[ind[1:sepPoint], 17:164])), "; partII(", ind[-(1:sepPoint)], ") median is", median(unlist(rawexp[ind[-(1:sepPoint)], 17:164])),
        "\nready for test!\n")
        
        testProbes <- ind[-(1:sepPoint)]
        #if(verbose) cat("I'm testpart, I have p", testProbes, "\n")
        restProbes <- c(ind[1:sepPoint], restProbes)
        #if(verbose) cat("I'm restpart, I have p", restProbes, "\n")
        
        #NOTE: min(ind[sepPoint], ind[sepPoint+1]) <- the probe just before the gap, for making plot.
        #      in 5'AS, it is the last probe of higher part in the first exon; in 3'AS, it is the last probe of the lower part in the last exon
        res <- ind[sepPoint]
        #>= P probes left in each group of the first exon, remember the gene name and do wilcox.test
        for(env in 1:4){
          ind_env <- which(as.numeric(menvironment) == env)
          res <- c(res, testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env))
          #if(verbose) cat("env", env, ":", testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env), "\n")
        }
        return(max(res))
      } else if(verbose) cat(" =>first exon increases/not enough exonProbes in each part/not median>=5 + median<5 after separation T^T\n")
    } else if(verbose) cat(" after permutation, no expExons T^T\n")
  } else if(verbose) cat(filename, "don't have enough exons/exonProbes/exonProbes in firstExon T^T\n")
}
#permutation5AS(filename, P=2, verbose=TRUE)


#get all expGenenames
allExpGenes <- NULL
for(chr in 1:5){
  allExpGenes <- c(allExpGenes, expGeneList[[chr]])
}


testRange=5000
maxList <- NULL
for(r in 1:1000){
  st <- proc.time()[3]
  
  #randomly select ** genes from all expGenes for each round of permutation
  genenames <- allExpGenes[sample(1:16448, testRange, replace=FALSE)]
  #filename = "AT1G01010"
  
  n=1
  while(n <= testRange){
    filename <- genenames[n]
    res <- permutation5AS(filename, P=2)
    if(!is.null(res)){
      maxList <- c(maxList, res)
      cat(filename, "is test in round", r, "-", n, "\n")
    }
    if(length(maxList) == 50*r) n=testRange+1
    else n <- n+1
  }
  if(length(maxList) < 50*r) maxList[(length(maxList)+1):(50*r)] <- NA
  
  et <- proc.time()[3]
  cat("permutation round", r, "is finished in", et-st, "s, and get values:", maxList[(50*r-49):(50*r)], "!\n\n")
}
write.table(maxList, file="Data/AS/Permutation for 5'AS.txt", row.names=FALSE)


png(filename="Data/AS/Histogram of Permutation for 5'AS.png", width=480, height=480, bg="White")
hist(maxList, breaks=50, border=grey(0.5), col="grey", main="Histogram of Results of Permutation for Alternative 5' Splicing Site", xlab="")
abline(v=59, col='red')
axis(1, at=59, labels=59, col='red', col.axis='red')
dev.off()
