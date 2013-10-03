#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 30-08-2013
# first written: 30-08-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#****************************************************** this is the final version for testing retained introns! ^_^ ******************************************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#
#main idea: compare each intron with all exons in this gene, no test whether intron expressed higher than 5(thre)!!!

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
testDffBtwParts <- function(exp_data=rawexp[ ,ind_env+16], testProbes, restProbes, verbose=FALSE){
  testPart <- as.numeric(unlist(exp_data[testProbes, ]))
  if(verbose) cat("We are testProbes:", testProbes, "\n")
  restPart <- as.numeric(unlist(exp_data[restProbes, ]))
  if(verbose) cat("We are restProbes:", restProbes, "\n")
  return(wilcox.test(testPart, restPart, alternative="less") $ p.value)# in RI, here's P-VALUE!!!
}
#testDffBtwParts(exp_data=rawexp[ ,ind_env+16], testProbes, restProbes, verbose=TRUE)

#intron retention test
permuteRI <- function(filename, P=2, verbose=FALSE){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueIntron <- unique(grep("intron", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have introns:", uniqueIntron, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")
  
  #for testing intron retention, at least 1 intron in a gene!!!
  if(length(uniqueIntron) > 0){
    if(verbose) cat(filename, "has", length(uniqueIntron), "introns, permute all the probes!\n")
    
    #shuffle all the labels in the gene, mix introns and exons
    randomOrder <- sample(probes_dir, length(probes_dir), replace=FALSE)
    rawexp[probes_dir, "tu"] <- rawexp[randomOrder, "tu"]
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
    if(verbose) cat("new tu order:", as.character(rawexp[probes_dir, "tu"]), "\n")
    
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
        res <- rep(NA, length(uniqueIntron))
      }else{
        if(verbose) cat("*in Env", env, "bgSet has", length(bg), "exonProbes, continue with test!\n")
        
        res <- NULL
        for(intron in uniqueIntron){
          inPrb <- probes_dir[rawexp[probes_dir, "tu"] == intron]
          
          if(length(inPrb) < P){
            if(verbose) cat("\t", intron, "has", length(inPrb), "intronProbes, <", P, "not enough probes to test T^T\n")
            res <- c(res, NA)
          }else{
            if(verbose) cat("\t", intron, "has", length(inPrb), "intronProbes, continue with wilcox.test!\n")
            res <- c(res, testDffBtwParts(exp_data=rawexp[ ,ind_env+16], testProbes=inPrb, restProbes=bg, verbose=FALSE))
          }
        }
      }
      resmatrix <- cbind(resmatrix, res)
    }
    return(resmatrix)
  }else if(verbose) cat(filename, "has no intron T^T\n")
}
#permuteRI(filename, P=2, verbose=FALSE)


#permute RI for chr 1-5
res <- NULL
for(n in 1:5){
  st <- proc.time()[3]
  cat("permutation round", n, "starts...\n")
  for(chr in 1:5){
    location <- paste0("Data/Raw/chr", chr, "_norm_hf_cor/")
    genenames <- gsub(".txt", "", dir(location)[grepl(".txt", dir(location)) & !grepl("_", dir(location))])
    
    pv = 0; nm = 1
    while(pv < 20){
      #filename = "AT1G01010"
      filename <- sample(genenames, 1, replace=FALSE)
      crntres <-  min(permuteRI(filename, P=2),na.rm=TRUE)
      if(is.finite(crntres)){
        res <- c(res, crntres)
        pv <- pv+1
        cat(filename, "\t")
      }else{
        cat("Selected a gene that doesn't allow for testing intron retention")
      }
      nm <- nm+1
    }
    cat("used", nm, "genes in total\n")
  }
  et <- proc.time()[3]
  cat("permutation round", n, "finished in", et-st, "s\n\n")
}
write.table(res, file=paste0("Data/AS/retainedIntron_permute.txt")) #***************************** change with goal !!! *****************************#

png(filename="Data/AS/Histogram of Permutation for RI.png", width=480, height=480, bg="White")
hist(res, breaks=50, border=grey(0.5), col="grey", main="Histogram of Results of Permutation for Retained Intron", xlab="", xlim=range(aa), ylim=c(0,65000))
abline(v=31.11, col='red')
axis(1, at=31.11, labels=31.11, col='red', col.axis='red')
dev.off()
