#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 19-06-2013
# first written: 21-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#********************************************** This is used to find alternative splicing at 5' or 3' **********************************************#


setwd("D:/Arabidopsis Arrays")

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

#skipping exon test
splicingTest <- function(filename, verbose = FALSE, ...){
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
    
    for(testExon in uniqueExon){
      ind <- exonID[rawexp[exonID, "tu"] == testExon]
      #if(verbose) cat(testExon, "has probes", ind, "\n")
      
      #at least 6 probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
      if(length(ind) >= 3){
        if(verbose) cat("\t***I'm", testExon, ", has", length(ind), "good probes, test me for skipping exon!\n")
        
        testProbes <- ind
        #if(verbose) cat("I'm testpart, I have p", testProbes, "\n")
        restProbes <- exonID[!exonID %in% ind]
        #if(verbose) cat("I'm restpart, I have p", restProbes, "\n")
        
        #>= 3 probes left in each group, remember the gene name and do t.test
        if(length(testProbes) >= 3 && length(restProbes) >= 3){
          if(verbose) cat(" =>separate at p", max(ind), ", have >= 3 good probes in each group, ready for test!\n")
          rownameList <- c(rownameList, filename)
          res <- max(ind)
          for(env in 1:4){
            ind_env <- which(as.numeric(menvironment) == env)
            res <- c(res, testDffBtwParts(exp_data = rawexp[ ,17:164], testProbes, restProbes, ind = ind_env, verbose, ...))
            #if(verbose) cat("env", env, ":", testDffBtwParts(exp_data = rawexp[ ,17:164], testProbes, restProbes, ind = ind_env, ...), "\n")
          }
          resmatrix <- rbind(resmatrix, res)
        } else if(verbose) cat(" =>after grouping don't have enough probes in each part T^T\n")
      } else if(verbose) cat("\t***I'm", testExon, ", not enough probes T^T\n")
    }
    rownames(resmatrix) <- rownameList
    return(resmatrix)
  } else if(verbose) cat("we don't have enough exons T^T\n")
}
#splicingTest(filename, useForTest = unlist, whichTest = wilcox.test, alternative = "less", verbose = TRUE)


#*************************************************************** test part ***************************************************************#
#Before start, mean, median or all individuals, change!!!
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- splicingTest(filename, useForTest = unlist, whichTest = wilcox.test, alternative = "less")
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      cat(filename, "is tested\n")
    }
  }
  colnames(resmatrix) <- c("sepProbe", "6H", "Dry_AR", "Dry_Fresh", "RP")
  
  write.table(resmatrix, file=paste0("Data/skippingExon/SE_chr", chr, "_wt_less.txt"), sep="\t") #********** change!!! **********#
  
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}



#*************************************************************** results analysis ***************************************************************#
#nTU
seMatrix <- NULL
for(chr in 1:5){
  seMatrix <- rbind(seMatrix, read.table(paste0("Data/skippingExon/SE_chr", chr, "_wt_less.txt"), row.names=NULL))
}
seThre <- -log10(0.05/nrow(seMatrix)/4)# = 6.31; 25508 exons were tested * 4 Env
rm(seMatrix)

resSEMatrixTU <- NULL
resSEMatrixGENE <- NULL
for(chr in 1:5){
  sechr <- read.table(paste0("Data/skippingExon/SE_chr", chr, "_wt_less.txt"), row.names=NULL)
  seGeneList <- list()
  
  resTU <- NULL
  resGENE <- NULL
  for(env in 1:4){
    resTU <- c(resTU, length(which(sechr[ ,env+2] >= seThre)))
    seGeneList[[paste0("env", env)]] <- unique(sechr[sechr[ ,env+2] >= seThre, 1])
    resGENE <- c(resGENE, length(seGeneList[[paste0("env", env)]]))
  }
  
  cnt_mixEnv <- 0
  tusGN <- NULL
  for(e in 1:nrow(sechr)){
    if(any(sechr[e, 3:6] >= seThre)){
      cnt_mixEnv <- cnt_mixEnv + 1
      tusGN <- c(tusGN, sechr[e,1])
    }
  }
  seGeneList$mixEnv <- unique(tusGN)
  
  resSEMatrixTU <- rbind(resSEMatrixTU, c(resTU, cnt_mixEnv))
  resSEMatrixGENE <- rbind(resSEMatrixGENE, c(resGENE, length(seGeneList$mixEnv)))
  save(seGeneList, file = paste0("Data/skippingExon/SE_chr", chr, "_genenameList.Rdata"))
}
rownames(resSEMatrixTU) <- paste0("chr", 1:5)
colnames(resSEMatrixTU) <- c(paste0("Env", 1:4), "mixEnv")
resSEMatrixTU

rownames(resSEMatrixGENE) <- paste0("chr", 1:5)
colnames(resSEMatrixGENE) <- c(paste0("Env", 1:4), "mixEnv")
resSEMatrixGENE
