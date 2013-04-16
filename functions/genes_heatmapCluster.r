#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 15-04-2013
# first written: 15-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]


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


#*************************************************************** load part ***************************************************************#
load(file="Data/fullModeMapping/expGenes.Rdata")


#************************************************************ make matrix part ************************************************************#
newprobematrix <- NULL
rownameList <- NULL
colnameList <- NULL

for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  genenames <- expGeneList[[chr]]
  
  #filename="AT1G01010.txt"
  for(filename in genenames){
    #cat(filename, "\n")
    rawexp <- read.table(paste0(location, filename), row.names=1, header=T)
    
    probes_dir <- probesDir(rawexp)
    exonID <- probes_dir[which(grepl("tu", rawexp[probes_dir, "tu"]))]
    
    newprobe <- NULL
    
    #Levels-numbers: 6H-35, Dry_AR-38, Dry_Fresh-39, RP-36
    for(env in 1:4){
      ind_env <- which(as.numeric(menvironment) == env)
      
      #make the colname of the right order for just one time
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
  rownameList <- c(rownameList, gsub(".txt", "", genenames))
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s!\n")
}
rownames(newprobematrix) <- rownameList
colnames(newprobematrix) <- colnameList
write.table(newprobematrix, file="Data/fullModeMapping/expGenes_1gene1probe.txt")

#to load this file
read.table("Data/fullModeMapping/expGenes_1gene1probe.txt", row.names=1, header=T)
