#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 08-07-2013
# first written: 08-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#******************************************************** this is the final version for testing AS at 5/3 site ^_^ *******************************************************#
#*********************************************************************** testing algorithm: ANOVA! ***********************************************************************#
#main idea:
#minimum of 6 probes in this exon
#Test how to split (using highest difference between two groups)
#Split into two groups
#test if every group has 3 probes
#YES -> T-Test the test group (5" start) against all the other probes in the gene
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
    dff <- sum(apply(exp_data[toGroup[length(toGroup):n], ], 2, median) - apply(exp_data[toGroup[1:(n-1)], ], 2, median))
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

#3'/5' AS test
#annotation: "goal" could be "5'AS" and "3'AS"
splicingTest35 <- function(filename, goal="5'AS", verbose=FALSE){
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
  if(length(uniqueExon) >= 2 && length(exonID) >= 7){
    if(verbose) cat(filename, "have", length(uniqueExon), "exons, and", length(exonID), "exonProbes!\n")
    
    if(goal == "5'AS") testExon <- uniqueExon[1]
    if(goal == "3'AS") testExon <- uniqueExon[length(uniqueExon)]
    
    ind <- exonID[rawexp[exonID, "tu"] == testExon]
    #if(verbose) cat(testExon, "has probes", ind, "\n")
    
    #at least 6 probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
    if(length(ind) >= 6){
      if(verbose) cat("\t***I'm", testExon, ", has", length(ind), "good probes, can be tested for", goal, "!\n")
      
      if(goal == "3'AS") ind <- sort(ind, decreasing=TRUE) #for last exon, overturn the order of probes!
      
      sepPoint <- findSepPoint(toGroup=ind, exp_data=rawexp[ ,17:164], verbose=FALSE)
      if(!is.null(sepPoint)){
        testProbes <- ind[1:sepPoint]
        #if(verbose) cat("I'm testpart, I have p", testProbes, "\n")
        restProbes <- exonID[!exonID %in% testProbes]
        #if(verbose) cat("I'm restpart, I have p", restProbes, "\n")
      } else{
        if(verbose) cat(" =>expression decreased at terminal\n")
        return()
      }
      
      #>= 3 probes left in each group of the first exon, remember the gene name and do t.test
      if(length(testProbes) >= 3 && length(ind[-(1:sepPoint)]) >= 3){
        if(verbose) cat(" =>separate between p", ind[sepPoint], "and p", ind[sepPoint+1], ", each group has >= 3 good probes, ready for test!\n")
        
        #NOTE: min(ind[sepPoint], ind[sepPoint+1]) <- the probe just before the gap, for making plot.
        #      in 5'AS, it is the last probe of higher part in the first exon; in3'AS, it is the last probe of the lower part in the last exon
        res <- min(ind[sepPoint], ind[sepPoint+1])
        for(env in 1:4){
          ind_env <- which(as.numeric(menvironment) == env)
          res <- c(res, testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env, verbose))
          #if(verbose) cat("env", env, ":", testDffBtwParts(exp_data = rawexp[ ,17:164], testProbes, restProbes, ind = ind_env), "\n")
        }
        return(res)
      } else if(verbose) cat(" =>after grouping don't have enough probes in each part T^T\n")
    } else if(verbose) cat("\t***I'm", testExon, ", not enough probes T^T\n")
  } else if(verbose) cat("we don't have enough exons T^T\n")
}
#splicingTest35(filename, goal="5'AS", verbose=TRUE)


#test 5'AS for chr 1-5
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  rownameList <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- splicingTest35(filename, goal="5'AS")
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      rownameList <- c(rownameList, filename)
      cat(filename, "is tested\n")
    }
  }
  if(!is.null(resmatrix)){
    rownames(resmatrix) <- rownameList
    colnames(resmatrix) <- c("sepProbe", "6H", "Dry_AR", "Dry_Fresh", "RP")
    write.table(resmatrix, file=paste0("Data/AS/5AS_chr", chr, "_ANOVA_less.txt"), sep="\t") #********** change!!! **********#
  } else cat("\tNO TEST!\n")
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}

#test 3'AS for chr 1-5
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  rownameList <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- splicingTest35(filename, goal="3'AS")
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      rownameList <- c(rownameList, filename)
      cat(filename, "is tested\n")
    }
  }
  if(!is.null(resmatrix)){
    rownames(resmatrix) <- rownameList
    colnames(resmatrix) <- c("sepProbe", "6H", "Dry_AR", "Dry_Fresh", "RP")
    write.table(resmatrix, file=paste0("Data/AS/3AS_chr", chr, "_ANOVA_less.txt"), sep="\t") #********** change!!! **********#
  } else cat("\tNO TEST!\n")
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}
