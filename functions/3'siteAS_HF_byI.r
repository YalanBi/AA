#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 27-08-2013
# first written: 27-08-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#********************************************* this is the final version for testing interaction regulated AS at 3 site ^_^ **********************************************#
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
#load genotype file
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)
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
findSepPoint <- function(toGroup=ind, exp_data=rawexp[ ,17:164], P=2, verbose=FALSE){
  dffs <- NULL
  for(n in (P+1):(length(toGroup)-P+1)){
    dff <- sum(apply(exp_data[toGroup[length(toGroup):n], ], 2, median)-apply(exp_data[toGroup[1:(n-1)], ], 2, median))
    dffs <- c(dffs, dff)
    if(verbose) cat("difference between probe(", toGroup[length(toGroup):n], ") and probe(", toGroup[1:(n-1)], ")is", dff, "\n")
  }
  if(min(dffs) < 0){
    if(verbose) cat("so dffList is:", dffs, "\n", "and edge probes are p", toGroup[which.min(dffs)+P-1], "and p", toGroup[which.min(dffs)+P], ", dff =", min(dffs), "\n")
    return(which.min(dffs)+P-1)
  }else return() #we want the testPart is lower than the other part of this exon, otherwise it is decay/decrease
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

#3'site IAS test
test3siteIAS <- function(filename, threTest=11.6, P=2, verbose=FALSE){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  int <- read.table(paste0("Data/FullModel/chr", chr, "_norm_hf_cor_FM/", filename, "_FM_Int.txt"), row.names=1, header=TRUE)
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")
  
  #for interaction regulated alternative splicing 3'site, at least 2 exons in a gene!!!
  if(length(uniqueExon) < 2){
    if(verbose) cat(filename, "has", length(uniqueExon), "exons, not enough T^T\n")
    return(c(0, 0, rep(0, 8)))
  }else{
    if(verbose) cat(filename, "has", length(uniqueExon), "exons!\n")
    
    ind <- sort(exonID[rawexp[exonID, "tu"] == uniqueExon[length(uniqueExon)]], decreasing=TRUE)
    if(length(ind) >= 2*P){
      if(verbose) cat("last exon", uniqueExon[length(uniqueExon)], "has probes", ind, ";")
      
      sepPoint <- findSepPoint(toGroup=ind, exp_data=rawexp[ ,17:164], P=P, verbose=FALSE)
    }else{
      if(verbose) cat("last exon", uniqueExon[length(uniqueExon)], "has", length(ind), "probes, not enough T^T\n")
      return(c(0, 0, rep(0, 8)))
    }
    
    if(is.null(sepPoint)){
      if(verbose) cat("but no right sep point T^T\n")
      return(c(0, 0, rep(0, 8)))
    }else{
      if(verbose) cat("and sep point is after p", ind[sepPoint+1], "!\n")
      
      partQTL <- int[ind[-(1:sepPoint)], ]
      if(any(apply(partQTL >= threTest, 2, sum) > 0)){
        m <- which.max(apply(partQTL >= threTest, 2, sum))
        
        #NOTE: min(ind[sepPoint], ind[sepPoint+1]) <- the probe just before the gap, for making plot.
        #      in 5'AS, it is the last probe of higher part in the first exon; in 3'AS, it is the last probe of the lower part in the last exon
        res <- c(ind[sepPoint+1], m)
        if(verbose) cat("and topMarker is", m, ", continue test!\n")
        
        for(env in 1:4){
          ind_env <- which(as.numeric(menvironment) == env)
          envGT1 <- ind_env[ind_env %in% which(geno[ ,m] == 1)]
          envGT2 <- ind_env[ind_env %in% which(geno[ ,m] == 2)]
          
          if(length(envGT1) > 0 && length(envGT2) > 0){
            for(gt in 1:2){
              gtInEnv <- ind_env[ind_env %in% which(geno[ ,m] == gt)]
              
              #check the part before sepPoint before get background set, if the median >= 5, continue!
              if(median(unlist(rawexp[ind[1:sepPoint], gtInEnv+16])) < 5){
                if(verbose) cat("*in Env", env, "first half median =", median(unlist(rawexp[ind[1:sepPoint], gtInEnv+16])), "< 5, too low T^T\n")
                res <- c(res, 0)
              }else{
                if(verbose) cat("*in Env", env, "first half median =", median(unlist(rawexp[ind[1:sepPoint], gtInEnv+16])), ">= 5, continue to get bgSet!\n")
                
                #get bgSet for each genotype in each Env; bgSet---the combination of TUs, the median of which is >= 5 of this genotype and in this Env
                bg <- NULL
                for(tu_bg in uniqueExon[-length(uniqueExon)]){
                  ind_bg <- exonID[rawexp[exonID, "tu"] == tu_bg]
                  if(length(ind_bg) > 0 && median(unlist(rawexp[ind_bg, gtInEnv+16])) >= 5){
                    if(verbose) cat("\tin Env", env, ": put", tu_bg, "into bgSet: median =", median(unlist(rawexp[ind_bg, gtInEnv+16])), ">= 5\n")
                    bg <- c(bg, ind_bg)
                  }
                }
                if(length(bg) < P){
                  if(verbose) cat("*in Env", env, "last half of bgSet has", length(bg), "exonProbes, <", P, " not enough to continue with test T^T\n")
                  res <- c(res, 0)
                }else{
                  bg <- c(bg, ind[1:sepPoint])
                  if(verbose) cat("*in Env", env, "bgSet is p", bg, ", continue with test!\n")
                  res <- c(res, testDffBtwParts(exp_data=rawexp[ ,gtInEnv+16], testProbes=ind[-(1:sepPoint)], restProbes=bg, verbose=FALSE))
                }
              }
            }
          }else{
            if(verbose) cat("in one gt, there's no RILs in Env", env, "T^T\n")
            res <- c(res, 0, 0)
          }
        }
      }else{
        res <- c(ind[sepPoint+1], 0, rep(0, 8))
        if(verbose) cat("but no topMarker T^T\n")
      }
      return(res)
    }
  }
}
#test3siteIAS(filename, threTest=11.6, P=2, verbose=TRUE)


#test 3'IAS for chr 1-5
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- test3siteIAS(filename, threTest=11.6, P=2)
    resmatrix <- rbind(resmatrix, res)
    cat(filename, "is tested\n")
  }
  rownames(resmatrix) <- genenames
  colnames(resmatrix) <- c("sepProbe", "topMarker", "6H/gt1", "6H/gt2", "Dry_AR/gt1", "Dry_AR/gt2", "Dry_Fresh/gt1", "Dry_Fresh/gt2", "RP/gt1", "RP/gt2")
  write.table(resmatrix, file=paste0("Data/geneticsAS/splicing3'siteByI_chr", chr, "_wt_p2.txt"), sep="\t")
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}
