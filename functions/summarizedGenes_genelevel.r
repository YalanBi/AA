#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 01-07-2013
# first written: 15-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#************************************************* this is the final version for summarizing probes into gene level! ^_^ *************************************************#
#******************* one probe for one gene, and the order of RILs is: first all individuals in Env1, next all individuals in Env2, then Env3 and Env4 *******************#
#****************************************************** this is used for  H E A T M A P  C L U S T E R I N G  ! ! ! ******************************************************#
setwd("D:/Arabidopsis Arrays")
menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]
#load exp genes
load(file="Data/ExpGenes/expGenes_final.Rdata")

#direction selection
probesDir <- function(exp_data = rawexp){
  if(unique(exp_data[,"strand"]) == "sense"){
    direction_id <- which(exp_data[, "direction"] == "reverse")
  }
  if(unique(exp_data[,"strand"]) == "complement"){
    direction_id <- which(exp_data[, "direction"] == "forward")
  }
  return(direction_id)
}

#summarize probes in one gene into one probe, each gene has one probe only, chr1-5 all in one file
newprobematrix <- NULL
rownameList <- NULL
colnameList <- NULL
for(chr in 1:5){
  location <- paste0("Data/Raw/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  genenames <- expGeneList[[chr]]
  
  #filename="AT1G01010"
  for(filename in genenames){
    rawexp <- read.table(paste0(location, filename, ".txt"), row.names=1, header=TRUE)
    probes_dir <- probesDir(rawexp)
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
    
    newprobe <- NULL
    #Levels-numbers: 6H-35, Dry_AR-38, Dry_Fresh-39, RP-36
    for(env in 1:4){
      ind_env <- which(as.numeric(menvironment) == env)
      
      #make the colname as the same order as the env for just one time
      if(length(colnameList) != 148){
        colnameList <- c(colnameList, colnames(rawexp[17:164])[ind_env])
        cat("colnames of env", env, "got!\n")
      }
      
      for(i in (ind_env + 16)){
        newprobe <- c(newprobe, median(rawexp[exonID, i]))
      }
      #cat("RILs in env", env, "finished\n")
    }
    newprobematrix <- rbind(newprobematrix, newprobe)
  }
  rownameList <- c(rownameList, genenames)
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s!\n")
}
rownames(newprobematrix) <- rownameList
colnames(newprobematrix) <- colnameList
write.table(newprobematrix, file="Data/summarizedGene/expGenes_1gene1probe.txt")

#to load this file
read.table("Data/summarizedGene/expGenes_1gene1probe.txt", row.names=1, header=TRUE)
