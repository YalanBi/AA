#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 18-06-2013
# first written: 18-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#****************************************** this is used to find those expressed genes from all the genes ******************************************#
#****************************** Here are 2 strategies. The first one is simpler, the second one has more requirements ******************************#


setwd("D:/Arabidopsis Arrays")

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


#*************************************************************** first/ old strategy ***************************************************************#
#main idea: find expressed genes, of which the median of all exon RILs is higher than threshold(expThre = 5)
expGeneTestSimple <- function(filename, use = median, expThre = 5, verbose = FALSE){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  probes_dir <- probesDir(rawexp)
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
  
  #there are some gene having only one probe and it is of wrong direction, which means there's no exon probe of right direction
  if(length(exonID) > 0){
    expDegree <- use(apply(rawexp[exonID, 17:164], 1, unlist))
    if(expDegree >= expThre){
      if(verbose) cat(filename, "exp degree is", expDegree, ", higher than", expThre, "!\n")
      return(TRUE)
    } else return(FALSE) #it's not expressed
  } else return(FALSE) #cause no exon probe left
  
}
#expGenes <- expGeneTestSimple(filename, use = median, expThre = 5)

#if necessary, save the results into a file
expGeneList <- vector("list", 5)
for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
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


#************************************************************** second/ new strategy ***************************************************************#
#main idea:
# goodGenes <- function(filename, N, P = -1, M = -1, Q = -1, S = -1){
#  has P probes per exons which are expressed at least above M
#  if(is there Q out of these exons that has S probes with a QTL){
#   goodGenes <- c(goodGenes , filename)
#  }
# has N or more exons
# end
expGeneTest <- function(filename, dirSelect = TRUE, exonSelect = TRUE, M = -1, use = median, P = -1, N = -1, verbose = FALSE){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  probeID <- 1:nrow(rawexp)
  
  if(dirSelect){
    if(unique(rawexp[,"strand"]) == "sense"){
      probeID <- which(rawexp[, "direction"] == "reverse")
    }
    if(unique(rawexp[,"strand"]) == "complement"){
      probeID <- which(rawexp[, "direction"] == "forward")
    }
    if(verbose)cat("after dirSelection, probes:", probeID, "\n")
  } else if(verbose)cat("no dirSelection\n")
  
  if(exonSelect){
    probeID <- probeID[grepl("tu", rawexp[probeID, "tu"])]
    if(verbose)cat("after exonSelect, probes:", probeID, "\n")
  } else if(verbose)cat("no exonSelect\n")
  
  if(M != -1){
    probeID <- probeID[apply(rawexp[probeID, 17:ncol(rawexp)], 1, use) >= M]
    if(verbose)cat("after expProbeTest, probes:", probeID, "\n")
  } else if(verbose)cat("no expTest\n")
  
  if(P != -1){
    index <- NULL
    for(tu in unique(rawexp[probeID, "tu"])){
      if(length(which(rawexp[probeID,"tu"] == tu)) >= P){
        index <- c(index, which(rawexp[probeID,"tu"] == tu))
      }
    }
    probeID <- probeID[index]
    if(verbose)cat("after expTuTest, probes:", probeID, ", nExpTu", length(unique(rawexp[probeID, "tu"])), "\n")
    if(length(probeID) == 0) return(FALSE);
  } else if(verbose)cat("no expTuTest\n")
  
  if(length(unique(rawexp[probeID, "tu"])) >= N) return(TRUE)
  else return(FALSE)# Not enough exons
}
#expGeneTest(filename, M = 4, P = 3, N = 2)

#if necessary, save the results into a file
expGeneList <- vector("list", 5)
for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  genenames <- gsub(".txt", "", dir(location)[grepl(".txt", dir(location)) & !grepl("_QTL.txt", dir(location)) & !grepl("_SNPin.txt", dir(location))])
  #filename "AT1G01010"
  for(filename in genenames){
    #cat(filename, "testing now\n")
    if(expGeneTest(filename, M = 4, P = 3, N = 2)){
      expGeneList[[chr]] <- c(expGeneList[[chr]], filename)
    }
  }
  et <- proc.time()[3]
  cat("find", length(expGeneList[[chr]]), "expressed genes on chr", chr, "and finished in", et-st, "s\n")
}
save(expGeneList,  file="Data/ExpGenes/expGenes_complicated.Rdata")
