#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 27-09-2013
# first written: 27-09-2013
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

#Find maximum difference between Part A and Part B of a single exon
findSepPoint <- function(toGroup=ind, exp_data=rawexp[ ,17:164], P=2, verbose=FALSE){
  dffs <- NULL
  if(length(toGroup) < P*2) return(NULL)
  for(n in (P+1):(length(toGroup)-(P-1))){
    dff <- abs(sum(apply(exp_data[toGroup[length(toGroup):n], ], 2, median) - apply(exp_data[toGroup[1:(n-1)], ], 2, median)))
    dffs <- c(dffs, dff)
    if(verbose) cat("difference between probe(", toGroup[length(toGroup):n], ") and probe(", toGroup[1:(n-1)], ")is", dff, "\n")
  }
  if(verbose) cat("so dffList is:", dffs, "\n", "and edge probes are p", toGroup[which.max(dffs)+P-1], "and p", toGroup[which.max(dffs)+P], ", dff =", max(dffs), "\n")
  return(which.max(dffs)+(P-1))
}
#findSepPoint(toGroup=ind, exp_data=rawexp[ ,17:164], P=2, verbose=TRUE)

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

#5'site AS test
test5siteAS <- function(filename, P=2, verbose=FALSE){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  if(rawexp[1,13]=="compliment") rawexp <- rawexp[nrow(rawexp):1, ]
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")
  
  #for alternative splicing 3 & 5'site, at least 2 exons in a gene!!!
  obj <- NULL; namez <- NULL
  if(length(uniqueExon) < 2){
    if(verbose) cat(filename, "has", length(uniqueExon), "exons, not enough T^T\n")
  }else{
    if(verbose) cat(filename, "has", length(uniqueExon), "exons!\n")
    for(exon in uniqueExon){
      ind <- exonID[rawexp[exonID, "tu"] == exon]
      pvals3 <- NULL; pvals5 <- NULL
      if(length(ind) >= 2*P){
        if(verbose) cat("Exon", exon, "has probes", ind, ";")
        sepPoint <- findSepPoint(toGroup=ind, exp_data=rawexp[ ,17:164], P=P, verbose=FALSE)
        if(is.null(sepPoint)){
          if(verbose) cat("but no right sep point T^T\n")
        }else{
          if(verbose) cat("and sep point is after p", ind[sepPoint], "!\n")
          toSplitAt <- ind[sepPoint]
          toContinueAt <- ind[sepPoint+1]
          for(env in 1:4){
            ind_env <- which(as.numeric(menvironment) == env)
            GroupA <- unlist(rawexp[ind[1]:toSplitAt,ind_env+16]); GroupB <- unlist(rawexp[toContinueAt:ind[length(ind)],ind_env+16])
            pvals3 <- c(pvals3, wilcox.test(GroupA, GroupB, alternative="less")$p.value)
            pvals5 <- c(pvals5, wilcox.test(GroupB, GroupA, alternative="less")$p.value)
          }
          if(exon == uniqueExon[1]) pvals3 <- c(1, 1, 1, 1)
          if(exon == uniqueExon[length(uniqueExon)]) pvals5 <- c(1, 1, 1, 1)
          obj <- rbind(obj, c(pvals3,pvals5))
          #cat(filename, exon, "3'", pvals3, "\n")
          #cat(filename, exon, "5'", pvals5, "\n")
          namez <- c(namez, paste0(filename, "_", exon))
          }
      }
    }
  }
  if(!is.null(obj)){
    colnames(obj) <- c(paste(c("6H", "Dry_AR", "Dry_Fresh", "RP"), c("3'"),sep="_"), paste(c("6H", "Dry_AR", "Dry_Fresh", "RP"), c("5'"),sep='_'))
    rownames(obj) <- namez
    return(obj)
  }
}
#test5siteAS(filename, P=2, verbose=TRUE)


#test 5'AS for chr 1-5
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- test5siteAS(filename, P=2)
    resmatrix <- rbind(resmatrix, res)
    cat(filename, "is tested\n")
  }
  #rownames(resmatrix) <- genenames
  #colnames(resmatrix) <- c("sepProbe", "6H", "Dry_AR", "Dry_Fresh", "RP")
  write.table(resmatrix, file=paste0("Data/AS/splicing3&5_chr", chr, "_wt_p2.txt"), sep="\t")
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}
