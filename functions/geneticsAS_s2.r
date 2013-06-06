#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 06-06-2013
# first written: 05-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#************************************************************** basic part **************************************************************#
setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]
#load genotype file
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)


#direction selection
probesDir <- function(exp_data, minExpression = -1){
  if(unique(exp_data[,"strand"]) == "sense"){
    direction_id <- which(exp_data[, "direction"] == "reverse")
  }
  if(unique(exp_data[,"strand"]) == "complement"){
    direction_id <- which(exp_data[, "direction"] == "forward")
  }
  if(minExpression == -1) return(direction_id)
  expressed <- which(apply(exp_data[,-c(1:16)], 1, median) > minExpression)
  return(expressed[direction_id %in% expressed])
}


#*************************************************************** load part ***************************************************************#
load(file="Data/fullModeMapping/expGenes.Rdata")


#function
# goodGenes <- function(filename, N, P = -1, M = -1, Q = -1, S = -1){
#  has P probes per exons which are expressed at least above M
#  if(is there Q out of these exons that has S probes with a QTL){
#   goodGenes <- c(goodGenes , filename)
#  }
# has N or more exons
# end


#************************************************************* test and count ************************************************************#
#count main eQTL(qtl >= 8.0 && int < 11.6)
countMainQTL <- function(chr = 1, threshold_qtl = 8.0, threshold_int = 11.6, cutoffratio = 0.6){
  location <- paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  genenames <- expGeneList[[chr]]
  QTLGene <- NULL
  for(filename in genenames){
    goodExons <- 0
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename), row.names=1, header=T)
    qtl <- read.table(paste0(location, gsub(".txt", "_FM_QTL.txt", filename)), row.names=1, header=T)
    
    probes_dir <- probesDir(rawexp, minExpression = 2)
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
    uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
    cnt_tu <- 1
    if(length(uniqueExon) >= 2){ #At least 2 (or more) exons
      for(cnt_tu in 1:length(uniqueExon)){
        #cat("The exon:", uniqueExon[cnt_tu],"\n")
        #cat("All exons:", rawexp[exonID, "tu"],"\n")
        tuID <- exonID[rawexp[exonID, "tu"] == uniqueExon[cnt_tu]] #tuID <- ID of exon probes of current tu name
        if(length(tuID) >= 3){
          if(any(apply(qtl[tuID,] > threshold_qtl, 2,sum) >= 3)){# At least a marker with 3 probes showing QTL
            cat("(", uniqueExon[cnt_tu]," ",length(tuID),") ")
            goodExons <- goodExons + 1
            if(goodExons >= 2){
              QTLGene <- c(QTLGene, filename)
              cat(filename, "at MARKER", m," IN tu",uniqueExon[cnt_tu],"has", max(apply(qtl[tuID,] > threshold_qtl, 2, sum)), " probes > 8 and finished in", et-st, "s!\n")
            }
            m <- which.max(apply(qtl[tuID,], 2,sum))
            et <- proc.time()[3]
            
          }
        }
      }
    }
  }
  save(QTLGene, file=paste0("Data/countQTL/mainQTL_chr", chr, ".Rdata"))
  QTLGene
}
QTLGene <- countMainQTL()

for(chr in 1:5){
  countMainQTL(chr, threshold_qtl = 8.0, threshold_int = 11.6, cutoffratio = 0.6))
}

for(filename in QTLGene){
  #load rawexp file, qtl file and int file
  rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename), row.names=1, header=T)
  #cat("rawexp loading succeed!\n")
  qtl <- read.table(paste0(location, gsub(".txt", "_FM_QTL.txt", filename)), row.names=1, header=T)
  #cat("qtl loading succeed!\n")
  int <- read.table(paste0(location, gsub(".txt", "_FM_Int.txt", filename)), row.names=1, header=T)
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  probes_dir <- probesDir(rawexp, minExpression=6)
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
  m <- which.max(apply(qtl, 2,sum))
  cat(filename," ",m,"\n")
  group1 <- which(geno[,m] == 1)
  group2 <- which(geno[,m] == 2)
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  probes_dir <- probesDir(rawexp)
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
  for(cnt_tu in 1:length(uniqueExon)){
    tuID <- exonID[rawexp[exonID, "tu"] == uniqueExon[cnt_tu]] #tuID <- ID of exon probes of current tu name
    if(length(tuID) >= 3){
      #if(any(apply(qtl[tuID,] > threshold_qtl, 2,sum) >= 3)){# At least a marker with 3 probes showing QTL
        #cat("QTL at TU:", qtl[tuID, m],"    <-> ");
        #cat("QTL NOT at TU:", qtl[exonID[-which(rawexp[exonID, "tu"] == uniqueExon[cnt_tu])], m],"\n");
        cat(uniqueExon[cnt_tu]," at marker",m, " ", t.test(qtl[tuID, m], qtl[exonID[-which(rawexp[exonID, "tu"] == uniqueExon[cnt_tu])], m])$p.value,"\n")
      #}
    }
  }
}

plot(c(0,nrow(rawexp)), c(0,10),t='n')
for(x in 1: nrow(rawexp)){
  points(rep(x,2), c(mean(unlist(rawexp[x,16+ which(geno[,4]==1)])),mean(unlist(rawexp[x,16+ which(geno[,4]==2)]))), col=c("black",'red'))
}


cnt_tu <- 1
tuID <- exonID[rawexp[exonID, "tu"] == uniqueExon[cnt_tu]] #tuID <- ID of exon probes of current tu name
if(length(tuID) >= 3){
  any(apply(qtl[tuID,] > threshold_qtl, 2,sum) >= 3)
}

#*************************************************************** load part ***************************************************************#
threshold_qtl =8; threshold_int =11.6; cutoffratio = 0.6
for(chr in 1:5){
  load(file = paste0("Data/countQTL/mainQTL_chr", chr, ".Rdata"))
  location <- paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  resMatrix <- NULL
  rownameList <- NULL
  
  genenames <- QTLGene #in order not to change QTLGene.Rdata
  for(filename in genenames){
    #cat(filename, "starts...\n")
    #load rawexp file, qtl file and int file
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename), row.names=1, header=T)
    #cat(" rawexp loading succeed!\n")
    qtl <- read.table(paste0(location, gsub(".txt", "_FM_QTL.txt", filename)), row.names=1, header=T)
    #cat(" qtl loading succeed!\n")
    int <- read.table(paste0(location, gsub(".txt", "_FM_Int.txt", filename)), row.names=1, header=T)
    #cat(" int loading succeed!\n")
    
    #get probes' and exons' ID
    probes_dir <- probesDir(rawexp)
    #cat(" probes of right direction:", probes_dir, "\n")
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
    #cat(" exons of right direction:", exonID, "\n")
    
    #uniqueExon <- all tu names of exon probes
    uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
    #cat(" tu names:", uniqueExon, "\n")
    
    for(tu in uniqueExon){
      tuID <- exonID[rawexp[exonID, "tu"] == uniqueExon[tu]]
      
      #could be there's one exon containing only one probe and the probe is of wrong direction. So have to check length(tuID) > 0.
      if(length(tuID) > 0){
        #cat(uniqueExon[cnt_tu], "has >= 1 probes, check my QTL and Int!\n")
        
        m <- 1; continue <- TRUE
        while(m %in% 1:716 && continue){
          #if there are at least 60% probes of this tu are qtl >= 8 and int < 11.6, then we know we can test this tu whether is cassette exon or not
          if(length(which(qtl[tuID, m] >= threshold_qtl && int[tuID, m] < threshold_int))/length(tuID) >= cutoffratio){
            cat(filename, uniqueExon[cnt_tu], "marker", m, "has main eQTL!!! and quit from testing!\n")
            continue <- FALSE
          } else m <- m + 1
        }
        
        #if continue is false, it means this exon has main eQTL. so move on to find the marker showing highest QTL!
        if(!continue){
          genoMarker <- which.max(colSums(qtl[tuID, ]))
          cat("I'm the marker", genoMarker, "with max QTL\n")
          
          #test 4 conds and 2 genotypes separately, test 8 times in total
          for(env in 1:4){
            ind_env <- which(as.numeric(menvironment) == env)
            
            res <- NULL
            for(gt in 1:2){
              testPart <- unlist(rawexp[tuID, ind_env[geno[ind_env,genoMarker]==gt]+16])
              restPart <- unlist(rawexp[exonID[-which(rawexp[exonID, "tu"] == uniqueExon[cnt_tu])], ind_env[geno[ind_env,genoMarker]==gt]+16])
              
              res <- c(res, -log10(wilcox.test(testPart, restPart)$p.value))
            }
            cat(filename, tu, env, res, "\n")
            
            resMatrix <- rbind(resMatrix, res)
            rownameList <- c(rownameList, paste0(gsub(".txt", "", filename), "_", tu, "_env", env))
          }
        }
      }
    }
  }
  
}
