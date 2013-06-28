##
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 28-06-2013
# first written: 21-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#******************************************************** this is the final version for testing AS at 5/3 site ^_^ *******************************************************#
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
  } else return()
}
#findSepPoint(toGroup=ind, exp_data=rawexp[ ,17:164], verbose=TRUE)

#test the difference between 2 groups
#annotation: useForTest = unlist -> use all individuals to do test, better than mean/median
#            whichTest <- wilcox.test/ t.test; one-side!!!
testDffBtwParts <- function(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env, useForTest=unlist, whichTest=wilcox.test, alternative="less", verbose=FALSE){
  testPart <- apply(exp_data[testProbes, ind], 2, useForTest)
  if(verbose) cat("We are testProbes:", testProbes, "\n")
  restPart <- apply(exp_data[restProbes, ind], 2, useForTest)
  if(verbose) cat("We are restProbes:", restProbes, "\n")
  return(-log10(whichTest(testPart, restPart, alternative) $ p.value))
}
#testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env, useForTest=unlist, whichTest=wilcox.test, alternative="less", verbose=TRUE)

#3'/5' AS test
#annotation: "goal" could be "5'AS" and "3'AS"
splicingTest35 <- function(filename, goal="5'AS", verbose=FALSE, ...){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")
  
  #for 5'/3' AS, at least 2 exons in a gene!!!
  if(length(uniqueExon) >= 2){
    if(verbose) cat(filename, "have", length(uniqueExon), "exons, >= 2!\n")
    
    resmatrix <- NULL
    rownameList <- NULL
    
    if(goal == "5'AS") testExon <- uniqueExon[1]
    if(goal == "3'AS") testExon <- uniqueExon[length(uniqueExon)]
    
    ind <- exonID[rawexp[exonID, "tu"] == testExon]
    #if(verbose) cat(testExon, "has probes", ind, "\n")
    
    #at least 6 probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
    if(length(ind) >= 6){
      if(verbose) cat("\t***I'm", testExon, ", has", length(ind), "good probes, can be tested for", goal, "!\n")
      
      if(goal == "3'AS") ind <- sort(ind, decreasing=TRUE)
      
      sepPoint <- findSepPoint(toGroup=ind, exp_data=rawexp[ ,17:164], verbose=FALSE)
      if(!is.null(sepPoint)){
        testProbes <- ind[1:sepPoint]
        #if(verbose) cat("I'm testpart, I have p", testProbes, "\n")
        restProbes <- exonID[!exonID %in% testProbes]
        #if(verbose) cat("I'm restpart, I have p", restProbes, "\n")
      } else return()
      
      #>= 3 probes left in each group, remember the gene name and do t.test
      if(length(testProbes) >= 3 && length(restProbes) >= 3){
        if(verbose) cat(" =>separate between p", ind[sepPoint], "and p", ind[sepPoint+1], ", each group has >= 3 good probes, ready for test!\n")
        rownameList <- c(rownameList, filename)
        
        #min(ind[sepPoint], ind[sepPoint+1]) <- the probe just before the gap, for making plot
        res <- min(ind[sepPoint], ind[sepPoint+1])
        for(env in 1:4){
          ind_env <- which(as.numeric(menvironment) == env)
          res <- c(res, testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes, restProbes, ind=ind_env, verbose, ...))
          #if(verbose) cat("env", env, ":", testDffBtwParts(exp_data = rawexp[ ,17:164], testProbes, restProbes, ind = ind_env, ...), "\n")
        }
        resmatrix <- rbind(resmatrix, res)
      } else if(verbose) cat(" =>after grouping don't have enough probes in each part T^T\n")
    } else if(verbose) cat("\t***I'm", testExon, ", not enough probes T^T\n")
    
    rownames(resmatrix) <- rownameList
    return(resmatrix)
  } else if(verbose) cat("we don't have enough exons T^T\n")
}
#splicingTest35(filename, goal="5'AS", useForTest=unlist, whichTest=wilcox.test, alternative="less", verbose=TRUE)


#test for chr 1-5
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- splicingTest35(filename, goal="5'AS", useForTest=unlist, whichTest=wilcox.test, alternative="less")
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      cat(filename, "is tested\n")
    }
  }
  colnames(resmatrix) <- c("sepProbe", "6H", "Dry_AR", "Dry_Fresh", "RP")
  
  write.table(resmatrix, file=paste0("Data/53terminalAS/5AS_chr", chr, "_wt_less.txt"), sep="\t") #********** change!!! **********#
  
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}


#******************************************************************* S T O P @ H E R E ! ! ! *************************************************************************#
#********************************************** This is used to find alternative splicing at 5' or 3' **********************************************#

#to find max difference between each probe in exp, for grouping
findMAXgap <- function(toGroup = ind, exp_data = rawexp[ ,17:164], minVSmax, verbose = FALSE){
  dff <- NULL
  for(n in 2:length(toGroup)){
    dff <- c(dff, sum(exp_data[toGroup[n], ] - exp_data[toGroup[n-1], ]))
    if(verbose) cat("difference between p", toGroup[n], "and p", toGroup[n-1], "is", sum(exp_data[toGroup[n], ] - exp_data[toGroup[n-1], ]), "\n")
  }
  if(verbose) cat("so dff is:", dff, "and largest gap is between p", toGroup[which(dff == minVSmax(dff))], "and p", toGroup[which(dff == minVSmax(dff))+1], "=", minVSmax(dff), "\n")
  
  return(which(dff == minVSmax(dff)))
}
#findMAXgap(toGroup = ind, exp_data = rawexp[ ,17:164], minVSmax, verbose = T)
#if it is the first exon, minVSmax = min; if it is the last exon, minVSmax = max

#to find the maximum difference between the two parts in first/last exon 
#annotation: useForTest = unlist -> use all individuals to do test, better than mean/median
#            whichTest <- wilcox.test/ t.test; one-side!!!
testDffBtwParts <- function(exp_data = rawexp[ ,17:164], testProbes, restProbes, ind = ind_env, useForTest = unlist, whichTest = wilcox.test, alternative = "less", verbose = FALSE){
  testPart <- apply(exp_data[testProbes, ind], 2, useForTest)
  if(verbose) cat("We are testProbes:", testProbes, "\n")
  restPart <- apply(exp_data[restProbes, ind], 2, useForTest)
  if(verbose) cat("We are restProbes:", restProbes, "\n")
  return(-log10(whichTest(testPart, restPart, alternative) $ p.value))
}
#testDffBtwParts(exp_data = rawexp[ ,17:164], testProbes, restProbes, ind = ind_env, useForTest = unlist, whichTest = wilcox.test, alternative = "less", verbose = TRUE)

#5'/3' AS test
splicingTest <- function(filename, goal = "5'AS", verbose = FALSE, ...){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")
  
  #for 5'/3' AS, at least 2 exons in a gene!!!
  if(length(uniqueExon) >= 2){
    if(verbose) cat(filename, "have", length(uniqueExon), "exons, >= 2!\n")
    
    resmatrix <- NULL
    rownameList <- NULL
    
    if(goal == "5'AS"){
      testExon <- uniqueExon[1]
      minVSmax <- min
    }
    if(goal == "3'AS"){
      testExon <- uniqueExon[length(uniqueExon)]
      minVSmax <- max
    }
    
    ind <- exonID[rawexp[exonID, "tu"] == testExon]
    #if(verbose) cat(testExon, "has probes", ind, "\n")
      
    #at least 6 probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
    if(length(ind) >= 6){
      if(verbose) cat("\t***I'm", testExon, ", has", length(ind), "good probes, test me for", goal, "!\n")
      
      #to find max difference between each probe in exp, to group them
      separatePoint <- findMAXgap(toGroup = ind, exp_data = rawexp[ ,17:164], minVSmax, verbose)
      
      testProbes <- ind[1:separatePoint]
      #if(verbose) cat("I'm testpart, I have p", testProbes, "\n")
      restProbes <- exonID[!exonID %in% ind[1:separatePoint]]
      #if(verbose) cat("I'm restpart, I have p", restProbes, "\n")
      
      #>= 3 probes left in each group, remember the gene name and do t.test
      if(length(testProbes) >= 3 && length(restProbes) >= 3){
        separateProbe <- ind[separatePoint]
        if(verbose) cat(" =>separate at p", separateProbe, ", have >= 3 good probes in each group, ready for test!\n")
        rownameList <- c(rownameList, filename)
        res <- separateProbe
        for(env in 1:4){
          ind_env <- which(as.numeric(menvironment) == env)
          res <- c(res, testDffBtwParts(exp_data = rawexp[ ,17:164], testProbes, restProbes, ind = ind_env, verbose, ...))
          #if(verbose) cat("env", env, ":", testDffBtwParts(exp_data = rawexp[ ,17:164], testProbes, restProbes, ind = ind_env, ...), "\n")
        }
        resmatrix <- rbind(resmatrix, res)
      } else if(verbose) cat(" =>after grouping don't have enough probes in each part T^T\n")
    } else if(verbose) cat("\t***I'm", testExon, ", not enough probes T^T\n")
  rownames(resmatrix) <- rownameList
  return(resmatrix)
  } else if(verbose) cat("we don't have enough exons T^T\n")
}
#splicingTest(filename, goal = "5'AS", useForTest = unlist, whichTest = wilcox.test, alternative = "less", verbose = TRUE)





#*************************************************************** test part ***************************************************************#
#goal = 5'AS
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- splicingTest(filename, goal = "5'AS", useForTest = unlist, whichTest = wilcox.test, alternative = "less")
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      cat(filename, "is tested\n")
    }
  }
  colnames(resmatrix) <- c("sepProbe", "6H", "Dry_AR", "Dry_Fresh", "RP")
  
  write.table(resmatrix, file=paste0("Data/53terminalAS/5'AS_chr", chr, "_wt_less.txt"), sep="\t") #********** change!!! **********#
  
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}
#goal = 3'AS
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- splicingTest(filename, goal = "3'AS", useForTest = unlist, whichTest = wilcox.test, alternative = "less")
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      cat(filename, "is tested\n")
    }
  }
  colnames(resmatrix) <- c("sepProbe", "6H", "Dry_AR", "Dry_Fresh", "RP")
  
  write.table(resmatrix, file=paste0("Data/53terminalAS/3'AS_chr", chr, "_wt_less.txt"), sep="\t") #********** change!!! **********#
  
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}


#*************************************************************** results analysis ***************************************************************#
#5'AS analysis
as5Matrix <- NULL
for(chr in 1:5){
  as5Matrix <- rbind(as5Matrix, read.table(paste0("Data/53terminalAS/5'AS_chr", chr, "_wt_less.txt"), row.names=NULL))
}
#Bonferroni correction
as5Thre <- -log10(0.05/nrow(as5Matrix)/4) #=5.11; 1622 first exons were tested * 4 Env
rm(as5Matrix)

res5Matrix <- NULL
for(chr in 1:5){
  as5chr <- read.table(paste0("Data/53terminalAS/5'AS_chr", chr, "_wt_less.txt"), row.names=NULL)
  as5ResList <- list()
  
  res <- NULL
  for(env in 1:4){
    res <- c(res, length(which(as5chr[ ,env+2] >= as5Thre)))
    as5ResList[[paste0("env", env)]] <- as5chr[as5chr[ ,env+2] >= as5Thre, 1]
  }
  
  cnt_mixEnv <- 0
  as5ResList$mixEnv <- NULL
  for(e in 1:nrow(as5chr)){
    if(any(as5chr[e, 3:6] >= as5Thre)){
      cnt_mixEnv <- cnt_mixEnv + 1
      as5ResList$mixEnv <- c(as5ResList$mixEnv, as5chr[e,1])
    }
  }
  
  res5Matrix <- rbind(res5Matrix, c(res, cnt_mixEnv))
  save(as5ResList, file = paste0("Data/53terminalAS/5'AS_chr", chr, "_genenameList.Rdata"))
}
rownames(res5Matrix) <- paste0("chr", 1:5)
colnames(res5Matrix) <- c(paste0("Env", 1:4), "mixEnv")
res5Matrix


#3'AS analysis
as3Matrix <- NULL
for(chr in 1:5){
  as3Matrix <- rbind(as3Matrix, read.table(paste0("Data/53terminalAS/3'AS_chr", chr, "_wt_less.txt"), row.names=NULL))
}
#Bonferroni correction
as3Thre <- -log10(0.05/nrow(as3Matrix)/4) #=5.09; 1523 last exons were tested * 4 Env
rm(as3Matrix)

res3Matrix <- NULL
for(chr in 1:5){
  as3chr <- read.table(paste0("Data/53terminalAS/3'AS_chr", chr, "_wt_less.txt"), row.names=NULL)
  as3ResList <- list()
  
  res <- NULL
  for(env in 1:4){
    res <- c(res, length(which(as3chr[ ,env+2] >= as3Thre)))
    as3ResList[[paste0("env", env)]] <- as3chr[as3chr[ ,env+2] >= as3Thre, 1]
  }
  
  cnt_mixEnv <- 0
  as3ResList$mixEnv <- NULL
  for(e in 1:nrow(as3chr)){
    if(any(as3chr[e, 3:6] >= as3Thre)){
      cnt_mixEnv <- cnt_mixEnv + 1
      as3ResList$mixEnv <- c(as3ResList$mixEnv, as3chr[e,1])
    }
  }
  
  res3Matrix <- rbind(res3Matrix, c(res, cnt_mixEnv))
  save(as3ResList, file = paste0("Data/53terminalAS/3'AS_chr", chr, "_genenameList.Rdata"))
}
rownames(res3Matrix) <- paste0("chr", 1:5)
colnames(res3Matrix) <- c(paste0("Env", 1:4), "mixEnv")
res3Matrix






#*************************************************************** main idea ***************************************************************#
rr <- NULL
for(x in 1:(nrow(rawexp)-1)){
  rr <- c(rr, sum(abs(rawexp[x,17:164]-rawexp[(x+1),17:164])))
}

#minimum of 4 probes
#Test how to split (using highest difference)
#Split into groups
#test if every group has 2 probes
#YES -> T-Test against the first group (5" start) all the other probes in the gene
#NO -> continue
