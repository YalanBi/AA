#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 05-09-2013
# first written: 05-09-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#******************************************************** this is the final version for testing skipping exons ^_^ *******************************************************#
#****************************************************** skipping exon: cassette exons + the first/last spliced exon ******************************************************#
#********************************************************** testing algorithm: Wilcox.test; cutoff: permutation **********************************************************#
#main idea: expGenes >= 2 exons
#           each exons >= 3 probes
#           testExon: meadian < 5; restExons: median >= 5; wilcox.test

setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]
#load exp genes
#load(file="Data/ExpGenes/expGenes_final.Rdata")

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
testDffBtwParts <- function(exp_data=rawexp[ ,ind_env+16], testProbes, restProbes, verbose=FALSE){
  testPart <- as.numeric(unlist(exp_data[testProbes, ]))
  if(verbose) cat("We are testProbes:", testProbes, "\n")
  restPart <- as.numeric(unlist(exp_data[restProbes, ]))
  if(verbose) cat("We are restProbes:", restProbes, "\n")
  return(-log10(wilcox.test(testPart, restPart, alternative="less") $ p.value))
}
#testDffBtwParts(exp_data=rawexp[ ,ind_env+16], testProbes, restProbes, verbose=TRUE)

#skipping exon test
permuteSE <- function(filename, P=3, verbose=FALSE){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")
  
  #for skipping exons, at least 2 exons in a gene!!!
  if(length(uniqueExon) < 2){
    if(verbose) cat(filename, "has", length(uniqueExon), "exons, not enough T^T\n")
    return(0)
  } else{
    if(verbose) cat(filename, "has", length(uniqueExon), "exons, permute exon probes!\n")
    
    randomOrder <- sample(exonID, length(exonID), replace=FALSE)
    rawexp[exonID, "tu"] <- rawexp[randomOrder, "tu"]
    if(verbose) cat("new tu order:", as.character(rawexp[exonID, "tu"]), "\n")
    
    resmatrix <- NULL
    for(env in 1:4){
      ind_env <- which(as.numeric(menvironment) == env)
      
      #get background set for each env; bgSet---the combination of TUs, the median of which is >= 5 in this Env
      bg <- NULL
      for(tu in uniqueExon){
        ind <- exonID[rawexp[exonID, "tu"] == tu]
        
        if(length(ind) > 0 && median(unlist(rawexp[ind, ind_env+16])) >= 5){
          if(verbose) cat("\tin Env", env, ": put", tu, "into bgSet: median =", median(unlist(rawexp[ind, ind_env+16])), ">= 5\n")
          bg <- c(bg, ind)
        }
      }
      
      if(length(bg) < P){
        if(verbose) cat("*in Env", env, "bgSet has", length(bg), "exonProbes, <", P, " not enough to continue with test T^T\n")
        res <- rep(0, length(uniqueExon))
      } else{
        if(verbose) cat("*in Env", env, "bgSet has", length(bg), "exonProbes, continue with test!\n")
        
        res <- NULL
        for(tu in uniqueExon){
          ind <- exonID[rawexp[exonID, "tu"] == tu]
          
          if(length(ind) < P){
            if(verbose) cat("\t", tu, "has", length(ind), "exonProbes, <", P, "not enough probes to test T^T\n")
            res <- c(res, 0)
          } else{
            if(verbose) cat("\t", tu, "has", length(ind), "exonProbes, continue with wilcox.test!\n")
            res <- c(res, testDffBtwParts(exp_data=rawexp[ ,ind_env+16], testProbes=ind, restProbes=bg, verbose=FALSE))
          }
        }
      }
      resmatrix <- cbind(resmatrix, res)
    }
    return(max(resmatrix))
  }
}
#permuteSE(filename, P=3, verbose=TRUE)


#permute SE for chr 1-5
res <- NULL
for(n in 1:1000){
  st <- proc.time()[3]
  cat("permutation round", n, "starts...\n")
  
  for(chr in 1:5){
    location <- paste0("Data/Raw/chr", chr, "_norm_hf_cor/")
    genenames <- gsub(".txt", "", dir(location)[grepl(".txt", dir(location)) & !grepl("_", dir(location))])
    testgenes <- sample(genenames, 20, replace=FALSE)
    cat(testgenes, "\n")
    
    #filename = "AT1G01010"
    for(filename in testgenes){
      res <- c(res, permuteSE(filename, P=3))
    }
  }
  et <- proc.time()[3]
  cat("permutation round", n, "finished in", et-st, "s\n\n")
}
write.table(res, file=paste0("Data/AS/skippingExon_permute.txt")) #***************************** change with goal !!! *****************************#

png(filename="Data/AS/Histogram of Permutation for SE.png", width=480, height=480, bg="White")
hist(res, breaks=50, border=grey(0.5), col="grey", main="Histogram of Results of Permutation for Skipping Exon", xlab="", xlim=range(aa), ylim=c(0,65000))
abline(v=31.11, col='red')
axis(1, at=31.11, labels=31.11, col='red', col.axis='red')
dev.off()
