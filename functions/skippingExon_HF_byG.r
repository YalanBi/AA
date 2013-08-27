#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 27-08-2013
# first written: 19-08-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#********************************************* this is the final version for testing skipping exons regulated by genetics^_^ *********************************************#
#****************************************************** skipping exon: cassette exons + the first/last spliced exon ******************************************************#
#******************************************************************** testing algorithm: Wilcox.test! ********************************************************************#
#main idea: expGenes >= 2 exons
#           each exons >= 3 probes
#           testExon: meadian < 5; restExons: median >= 5; wilcox.test

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

#G/I regulated skipping exon test
testSEbyG <- function(filename, threTest=8, P=3, verbose=FALSE){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  qtl <- read.table(paste0("Data/FullModel/chr", chr, "_norm_hf_cor_FM/", filename, "_FM_QTL.txt"), row.names=1, header=TRUE)
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")
  
  #for skipping exons by G/I, at least 2 exons in a gene!!!
  if(length(uniqueExon) < 2){
    if(verbose) cat(filename, "has", length(uniqueExon), "exons, not enough T^T\n")
    return(matrix(c(uniqueExon, 0, rep(0, 8)), ncol=10, dimnames=list(filename)))
  }else{
    if(verbose) cat(filename, "has", length(uniqueExon), "exons!\n")
    
    resmatrix <- NULL
    for(tu in uniqueExon){
      ind <- exonID[rawexp[exonID, "tu"] == tu]
      tuQTL <- qtl[ind, ]
      
      if(length(ind) >= P && any(apply(tuQTL >= threTest, 2, sum) > 0)){
        m <- which.max(apply(tuQTL >= threTest, 2, sum))
        res <- c(tu, m)
        if(verbose) cat("", tu, "has topMarker", m, ", continue test!\n")
        
        for(env in 1:4){
          ind_env <- which(as.numeric(menvironment) == env)
          envGT1 <- ind_env[ind_env %in% which(geno[ ,m] == 1)]
          envGT2 <- ind_env[ind_env %in% which(geno[ ,m] == 2)]
          
          if(length(envGT1) > 0 && length(envGT2) > 0){
            for(gt in 1:2){
              gtInEnv <- ind_env[ind_env %in% which(geno[ ,m] == gt)]
              
              #get bgSet for each genotype in each Env; bgSet---the combination of TUs, the median of which is >= 5 of this genotype and in this Env
              bg <- NULL
              for(tu_bg in uniqueExon){
                ind_bg <- exonID[rawexp[exonID, "tu"] == tu_bg]
                if(length(ind_bg) > 0 && median(unlist(rawexp[ind_bg, gtInEnv+16])) >= 5){
                  if(verbose) cat("\tin Env", env, "+ gt", gt, ": put", tu_bg, "into bgSet: median =", median(unlist(rawexp[ind_bg, gtInEnv+16])), ">= 5\n")
                  bg <- c(bg, ind_bg)
                }
              }
              if(length(bg) < P){
                if(verbose) cat("*in Env", env, "bgSet for gt", gt, "has", length(bg), "exonProbes, <", P, " not enough to continue with test T^T\n")
                res <- c(res, 0)
              }else{
                if(verbose) cat("*in Env", env, "bgSet for gt", gt, "has", length(bg), "exonProbes, continue with test!\n")
                res <- c(res, testDffBtwParts(exp_data=as.matrix(rawexp[ ,gtInEnv+16]), testProbes=ind, restProbes=bg, verbose=FALSE))
              }
            }
          }else{
            if(verbose) cat("in one gt, there's no RILs in Env", env, "T^T\n")
            res <- c(res, 0, 0)
          }
        }
        resmatrix <- rbind(resmatrix, res)
      }else{
        if(verbose) cat("", tu, "has no topMarkers T^T\n")
        resmatrix <- rbind(resmatrix, c(tu, 0, rep(0, 8)))
      }
    }
    rownames(resmatrix) <- rep(filename, length(uniqueExon))
    return(resmatrix)
  }
}
#testSEbyG(filename, threTest=8, P=3, verbose=TRUE)

for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts testing eQTL regulated skipping exons...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    
    resmatrix <- rbind(resmatrix, testSEbyG(filename, threTest=8, P=3))
    cat(filename, "is tested\n")
  }
  colnames(resmatrix) <- c("sepProbe", "topMarker", "6H/gt1", "6H/gt2", "Dry_AR/gt1", "Dry_AR/gt2", "Dry_Fresh/gt1", "Dry_Fresh/gt2", "RP/gt1", "RP/gt2")
  write.table(resmatrix, file=paste0("Data/geneticsAS/skippingExonByG_chr", chr, "_wt_p3.txt"), sep="\t")
  
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}



#qtl <- read.table(paste0("Data/summarizedGene/expGenes_chr", chr, "_FMD_QTL.txt"), row.names=1, header=TRUE)
#fnPart <- abs(qtl[grep(filename, rownames(qtl)), ])
