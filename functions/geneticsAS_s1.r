#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 31-05-2013
# first written: 31-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]


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


ce_threshold = 5.86; threshold_qtl = 8.0; threshold_int = 11.6; cutoffnProbe = 2; #cutoffratio = 0.6#(cutoffratio>=0.6; cutoffnProbe >= 2 probes)
genesGAS <- list()
genesIAS <- list()
for(chr in 1:5){
  location <- paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  cematrix <- read.table(paste0("Data/cassetteExon/cassetteExon_chr", chr, "_allind.txt"), row.names=1, header=T)
  
  #for all genes which have cassette exons in at least one env
  genenames <- sort(unique(unlist(lapply(strsplit(rownames(which(cematrix >= ce_threshold, arr.ind=T)), "_"), "[[", 1))))
  
  genesGAS[[paste0("chr", chr)]] <- vector("list",4)
  genesIAS[[paste0("chr", chr)]] <- vector("list", 4)
  
  #filename(AT1G01010)
  for(filename in genenames){
    #cat(filename, "starts!\n")
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
    #cat(" rawexp loading succeed!\n")
    qtl <- read.table(paste0(location, filename, "_FM_QTL.txt"), row.names=1, header=T)
    #cat(" qtl loading succeed!\n")
    int <- read.table(paste0(location, filename, "_FM_Int.txt"), row.names=1, header=T)
    #cat(" int loading succeed!\n")
    
    #get probes' and exons' ID
    probes_dir <- probesDir(rawexp)
    #cat("probes of right direction:", probes_dir, "\n")
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
    #cat("exons of right direction:", exonID, "\n")
    
    #uniqueExon <- all tu names of exon probes
    uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
    #cat("tu names:", as.character(uniqueExon), "\n")
    
    for(env in 1:4){
      CEinEnv <- unlist(lapply(strsplit(grep(filename, rownames(cematrix), value=TRUE)[cematrix[grep(filename, rownames(cematrix), value=TRUE),env]>=ce_threshold], "_"), "[[", 2))
      #cat("=>Env", env, "has cassette exons:", CEinEnv, "\n")
      
      for(ce in CEinEnv){
        tuID <- exonID[rawexp[exonID, "tu"] == ce]
        #cat(" *", ce, "have right direction probes", tuID, "\n")
        
        consSigQTLMarkers <- NULL
        consSigIntMarkers <- NULL
        for(m in 1:716){
          #************************************************* make sure the criterion ! ! ! *************************************************#
          #if(length(which(qtl[tuID, m] >= threshold_qtl && int[tuID, m] < threshold_int))/length(tuID) >= cutoffratio){
          if(length(which(qtl[tuID, m] >= threshold_qtl && int[tuID, m] < threshold_int)) >= cutoffnProbe){
            consSigQTLMarkers <- c(consSigQTLMarkers, m)
            #cat("  ^^marker", m, "has main eQTL in", tuID, "\n")
          }
          #if(length(which(int[tuID, m] >= threshold_int))/length(tuID) >= cutoffratio){
          if(length(which(int[tuID, m] >= threshold_int)) >= cutoffnProbe){
            consSigIntMarkers <- c(consSigIntMarkers, m)
            #cat("  ^^marker", m, "has Int QTL in", tuID, "\n")
          }
        }
        
        if(length(consSigQTLMarkers) > 0){
          genesGAS[[paste0("chr", chr)]][[env]] <- c(genesGAS[[paste0("chr", chr)]][[env]], paste0(filename, "_", ce))
          cat("  in env", env, ",", filename, "'s", ce, "has main eQTL, remember me!!!\n")
        } #else cat("  in env", env, ",", filename, "'s", ce, "has noooooooooo main eQTL T^T\n")
        if(length(consSigIntMarkers) > 0){
          genesIAS[[paste0("chr", chr)]][[env]] <- c(genesIAS[[paste0("chr", chr)]][[env]], paste0(filename, "_", ce))
          cat("  in env", env, ",", filename, "'s", ce, "has Int QTL, remember me!!!\n")
        } #else cat("  in env", env, ",", filename, "'s", ce, "has noooooooooo Int QTL T^T\n")
      }
    }
    
  }
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}
genesGAS
genesIAS
