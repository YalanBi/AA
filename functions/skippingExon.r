#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 02-05-2013
# first written: 02-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
#load environment file
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

#3'/5' skipping exon test
test5Skipping <- function(expdata = rawexp[exonID, ind_env + 16], ind, use){
  fexonProbe <- apply(expdata[ind, ], 2, use)
  #cat(exonprobe, "\n")
  otherProbe <- apply(expdata[!ind, ], 2, use)
  #cat(otherprobe, "\n")
  return(-log10(t.test(exonprobe, otherprobe, alternative="less")$p.value))
}


#*************************************************************** load part ***************************************************************#
#load exp genes
load(file="Data/fullModeMapping/expGenes.Rdata")


#*************************************************************** test part ***************************************************************#
cnt_mem <- NULL
#mean or median, change!!!
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  rownameList <- NULL
  
  cnt <- 0
  
  for(filename in genenames){
    cat(filename, "...\n")
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename), row.names=1, header=T)
    #cat("rawexp loading succeed!\n")
    
    #get probes' and exons' ID
    probes_dir <- probesDir(rawexp)
    #cat("probes of right direction:", probes_dir, "\n")
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
    #cat("exons of right direction:", exonID, "\n")
    
    #uniqueExon <- all tu names of exon probes
    uniqueExon <- unique(rawexp[exonID, "tu"])
    #cat("tu names:", as.character(uniqueExon), "\n")
    
    #at least 2 exons in a gene!!!
    if(length(uniqueExon) >= 2){
      #cat("we have >= 2 exons!\n")
      gene_count <- gene_count+1
      exon <- uniqueExon[1]
        #ind <- judge which probe in exonID is of current exon name (T/F)
        ind <- rawexp[exonID, "tu"] == exon
        #cat(as.character(exon), "has probes", exonID[ind], "\n")
        if(length(which(ind)) >= 2){
          #cat("I'm", exon, ">= 2 good probes, t.test me and remember my tu name!\n")
          rownameList <- c(rownameList, paste0(gsub(".txt", "", filename), "_", exon))
          
          res <- NULL
          for(env in 1:4){
            ind_env <- which(as.numeric(menvironment) == env)
            
            #use mean/median to test cassette
            res <- c(res, testCassette(expdata = rawexp[exonID, ind_env + 16], ind, use = mean)) #********** change!!! **********#
            
            #t.test once, counter plus 1!!! to calculate the 
            cnt <- cnt + 1
          }
       resmatrix <- rbind(resmatrix, res)
        }

    }
  }
  cnt_mem <- c(cnt_mem, cnt)
  
  rownames(resmatrix) <- rownameList
  colnames(resmatrix) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
  write.table(resmatrix, file=paste0("Data/cassetteExon/cassetteExon_chr", chr, "_mean_D2p.txt"), sep="\t") #change!!! 
  
  et <- proc.time()[3]
  cat("and finished in", et-st, "s\n")
}