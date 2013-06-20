#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 18-06-2013
# first written: 18-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#******************************************** this is the final version, based on results from this --- ********************************************#
#****************************************************** --- continue all following analysis!! ******************************************************#
#main idea: find expressed genes, of which the median of all exon RILs is higher than threshold(expThre = 5) in any environment!

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

expGeneTestSimple <- function(filename, use = median, expThre = 5, verbose = FALSE){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  probes_dir <- probesDir(rawexp)
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
  
  #there are some gene having only one probe and it is of wrong direction, which means there's no exon probe of right direction
  if(length(exonID) > 0){
    #continue <- TRUE
    for(env in 1:4){
      ind_env <- which(as.numeric(menvironment) == env)
      expDegree <- use(apply(rawexp[exonID, ind_env+16], 1, unlist))
      if(expDegree >= expThre){
        if(verbose) cat(filename, "exp degree is", expDegree, ", higher than", expThre, "in env", env, "!\n")
        return(TRUE)
      } else if(verbose) cat("not expressed in env", env, "\n")
    }
    if(env == 4) return(FALSE)
  } else return(FALSE) #cause no exon probe left
}
#expGenes <- expGeneTestSimple(filename, use = median, expThre = 5, verbose = TRUE)

#if necessary, save the results into a file
expGeneList <- vector("list", 5)
for(chr in 1:5){
  location <- paste0("Data/Raw/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  genenames <- gsub(".txt", "", dir(location)[grepl(".txt", dir(location)) & !grepl("_QTL.txt", dir(location)) & !grepl("_SNPin.txt", dir(location))])
  #filename "AT1G01010"
  for(filename in genenames){
    cat(filename, "testing now\n")
    if(expGeneTestSimple(filename, use = median, expThre = 5)){
      expGeneList[[chr]] <- c(expGeneList[[chr]], filename)
    }
  }
  et <- proc.time()[3]
  cat("find", length(expGeneList[[chr]]), "expressed genes on chr", chr, "and finished in", et-st, "s\n")
}
save(expGeneList,  file="Data/ExpGenes/expGenes_simple.Rdata")



#check overlap
setwd("D:/Arabidopsis Arrays")
#load exp genes
load(file="Data/ExpGenes/expGenes_simple.Rdata")
sim <- expGeneList
rm(expGeneList)
load(file="Data/ExpGenes/expGenes_complicated_m5.Rdata")
com <- expGeneList
rm(expGeneList)
resm <- NULL
for(chr in 1:5){
  cat(length(sim[[chr]]), length(com[[chr]]), "\n")
  cat(length(which(sim[[chr]] %in% com[[chr]])), length(which(com[[chr]] %in% sim[[chr]])), "\n")
  cat(length(which(!sim[[chr]] %in% com[[chr]])), length(which(!com[[chr]] %in% sim[[chr]])), "\n\n")
  resm <- rbind(resm, c(length(sim[[chr]]), length(which(sim[[chr]] %in% com[[chr]])), length(com[[chr]])))
}
rownames(resm) <- 1:5
colnames(resm) <- c("old", "overlap", "new")
resm
