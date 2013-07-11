#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 11-07-2013
# first written: 08-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#****************************************************** this is the final version for testing retained introns! ^_^ ******************************************************#
#***************************************************************** the same algorithm as skipping exons! *****************************************************************#
#*********************************************************************** testing algorithm: ANOVA! ***********************************************************************#
#main idea: compare each intron with all exons in this gene, no test whether intron expressed higher than 5(thre)!!!

setwd("D:/Arabidopsis Arrays")
#load environment file
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

#intron retention test
splicingTestRI <- function(filename, P=2, verbose=FALSE){
  if(verbose) cat("now is testing retained introns!\n")
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueIntron <- unique(grep("intron", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have introns:", uniqueIntron, "\n")
  
  #for testing intron retention, at least P exon probes and 1 intron in a gene!!!
  if(length(exonID) >= P && length(uniqueIntron) >= 1){
    if(verbose) cat(filename, "have", length(exonID), "exon probes of right direction and", length(uniqueIntron), "introns!\n")
    
    resmatrix <- NULL
    rownameList <- NULL
    
    for(testIntron in uniqueIntron){
      ind <- probes_dir[rawexp[probes_dir, "tu"] == testIntron]
      #if(verbose) cat(testIntron, "has probes", ind, "\n")
      
      #at least P probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
      if(length(ind) >= P){
        if(verbose) cat("\t***I'm", testIntron, ", has", length(ind), "intron probes of right direction, test me for retained intron!\n")
        
        #>= P probes left in each group, remember the gene name and do wilcox.test
        rownameList <- c(rownameList, filename)
        res <- max(ind)
        for(env in 1:4){
          ind_env <- which(as.numeric(menvironment) == env)
          res <- c(res, testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes=ind, restProbes=exonID, ind=ind_env, verbose))
          #if(verbose) cat("env", env, ":", testDffBtwParts(exp_data=rawexp[ ,17:164], testProbes=ind, restProbes=exonID, ind=ind_env), "\n")
        }
        resmatrix <- rbind(resmatrix, res)
      } else if(verbose) cat("\t***I'm", testIntron, ", not enough intron probes T^T\n")
    }
    rownames(resmatrix) <- rownameList
    return(resmatrix)
  } else if(verbose) cat(filename, "we don't have any introns or not enough exon probes T^T\n")
}
#splicingTestRI(filename, P=2, verbose=TRUE)


#test for chr 1-5
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- splicingTestRI(filename, P=2)
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      cat(filename, "is tested\n")
    }
  }
  if(!is.null(resmatrix)){
    colnames(resmatrix) <- c("sepProbe", "6H", "Dry_AR", "Dry_Fresh", "RP")
    write.table(resmatrix, file=paste0("Data/AS/RI_chr", chr, "_ANOVA_less.txt"), sep="\t") #**************************** change with goal!!! ****************************#
  } else cat("\tNO TEST!\n")
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}
