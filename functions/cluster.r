#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 17-07-2013
# first written: 17-07-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#***************************************************** this is the final version for getting data for clustering ^_^ *****************************************************#
#****************************************************************** testing algorithm: t.test$statistic ******************************************************************#

setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]
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

#t.test each env against the other three
getClusterData <- function(filename, P=3, verbose=FALSE){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=TRUE)
  newexp <- rawexp[ , 17:164]
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  
  if(length(exonID) >= P){
    if(verbose) cat("I'm", filename, ", have more than", P, "probes\n")
    
    res <- NULL
    for(env in 1:4){
      ind_env <- which(as.numeric(menvironment) == env)
      res <- c(res, t.test(newexp[exonID, ind_env], newexp[exonID, -ind_env])$statistic)
      if(verbose) cat("env", env, "vs the other three, t is", t.test(newexp[exonID, ind_env], newexp[exonID, -ind_env])$statistic, "\n")
    }
    return(res)
  } else if(verbose) cat("I'm", filename, "less than", P, "probes...\n")
}


resmatrix <- NULL
rownameList <- NULL
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  for(filename in genenames){
    res <- getClusterData(filename, P=3)
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      rownameList <- c(rownameList, filename)
    } else cat(filename, "skipping, too few probes...\n")
  }
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s!\n")
}
rownames(resmatrix) <- rownameList
colnames(resmatrix) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
write.table(resmatrix, file="Data/cluster/dataForCluster.txt", sep="\t")



#get top 100 genes in each env, both high and low
resmatrix <- read.table("Data/cluster/dataForCluster.txt", row.names=1, header=TRUE)

topHigh <- NULL
for(env in 1:4){
  res <- sort(resmatrix[ ,env], decreasing=TRUE)
  topHigh <- cbind(topHigh, rownames(as.matrix(res[1:100])))
}
colnames(topHigh) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
write.table(topHigh, file="Data/cluster/top100+.csv", sep=" ", row.names=FALSE, col.names=TRUE)

topLow <- NULL
for(env in 1:4){
  res <- sort(resmatrix[ ,env], decreasing=FALSE)
  topLow <- cbind(topLow, rownames(as.matrix(res[1:100])))
}
colnames(topLow) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
write.table(topLow, file="Data/cluster/top100-.csv", sep=" ", row.names=FALSE, col.names=TRUE)
