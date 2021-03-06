#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 06-06-2013
# first written: 05-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#************************************************************** basic part **************************************************************#
setwd("D:/Arabidopsis Arrays")

#**************************************************** fisrt select potential QTLGene ****************************************************#
load(file="Data/fullModeMapping/expGenes.Rdata")

#function
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



potentialQTLGene <- function(filename, N = -1, toTest = "QTL", threTest = 8, Q = -1, S = -1, verbose = FALSE){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  
  if(unique(rawexp[,"strand"]) == "sense"){
    probeID <- which(rawexp[, "direction"] == "reverse")
  }
  if(unique(rawexp[,"strand"]) == "complement"){
    probeID <- which(rawexp[, "direction"] == "forward")
  }
  if(verbose)cat("after dirSelection, probes:", probeID, "\n")

  probeID <- probeID[grepl("tu", rawexp[probeID, "tu"])]
  if(verbose)cat("after exonSelect, probes:", probeID, "\n")
  
  hasQTL <- 0;
  if(length(unique(rawexp[probeID, "tu"])) >= N){
    if(verbose)cat("Starting QTL test with", length(unique(rawexp[probeID, "tu"])), "Exons\n");
    if(S != -1){
      testQTL <- read.table(paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/", filename, "_FM_", toTest, ".txt"), row.names=1, header=T)
      index <- NULL
      for(tu in unique(rawexp[probeID, "tu"])){
        tuID <- probeID[rawexp[probeID, "tu"] == tu]
        if(any(apply(testQTL[tuID,] > threTest, 2, sum) >= S)){
          index <- c(index, tuID)
          hasQTL <- hasQTL + 1;
        }
      }
      if(verbose)cat("Found", hasQTL, "Exons with at least", S, "probes on there (", index, ")\n");
      return(hasQTL >= Q)
    }else{ # Has enough Exons and no QTL check so its Good
      return(TRUE); #We dont care if S == -1
    }
  }else{ # Not enough exons
    return(FALSE)
  }
}


#************************************************************** t.test part **************************************************************#
#direction selection
probesDir <- function(rawexp, minExpression = -1){
  if(unique(rawexp[,"strand"]) == "sense"){
    direction_id <- which(rawexp[, "direction"] == "reverse")
  }
  if(unique(rawexp[,"strand"]) == "complement"){
    direction_id <- which(rawexp[, "direction"] == "forward")
  }
  
  if(minExpression != -1){
    expressed <- which(apply(rawexp[ ,17:164], 1, median) >= minExpression)
    return(direction_id[direction_id %in% expressed])
  }else return(direction_id)
}


#ttest exp of one exon against the rest, env and geno separately
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)

wilcoxtestQTLGene <- function(filename, toTest = "QTL", P = 3, verbose = FALSE){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  testQTL <- read.table(paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/", filename, "_FM_", toTest, ".txt"), row.names=1, header=T)
  
  probes_dir <- probesDir(rawexp, minExpression = -1)
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  
  resMatrix <- NULL
  resRN <- NULL
  for(tu in uniqueExon){
    tuID <- exonID[rawexp[exonID, "tu"] == tu] #tuID <- ID of exon probes of current tu name
    if(length(tuID) >= P){# && length(exonID[-which(rawexp[exonID, "tu"] == tu)]) >= P
      m <- which.max(apply(testQTL[tuID], 2, sum))
      if(verbose) cat(filename, tu, "most sig marker is ", m, "\n")
      geno1 <- which(geno[,m] == 1)
      geno2 <- which(geno[,m] == 2)
      
      res1vr <- NULL
      res2vr <- NULL
      #res1v2 <- NULL
      for(env in 1:4){
        ind_env <- which(as.numeric(menvironment) == env)
        envGroup1 <- ind_env[ind_env %in% geno1]
        envGroup2 <- ind_env[ind_env %in% geno2]
        if(length(envGroup1) != 0 && length(envGroup2) != 0){
          res1vr <- c(res1vr, -log10(wilcox.test(unlist(rawexp[tuID, envGroup1 + 16]), unlist(rawexp[exonID[-which(rawexp[exonID, "tu"] == tu)], envGroup1 + 16]), alternative = "less")$p.value))
          res2vr <- c(res2vr, -log10(wilcox.test(unlist(rawexp[tuID, envGroup2 + 16]), unlist(rawexp[exonID[-which(rawexp[exonID, "tu"] == tu)], envGroup2 + 16]), alternative = "less")$p.value))
          #res1v2 <- c(res1v2, -log10(wilcox.test(rawexp[tuID, envGroup1 + 16], rawexp[tuID, envGroup2 + 16])$p.value))
        }
        if(length(envGroup1) == 0){
          res1vr <- c(res1vr, -1)
          #res1v2 <- c(res1v2, -1)
        }
        if(length(envGroup2) == 0){
          res2vr <- c(res2vr, -2)
          #res1v2 <- c(res1v2, -2)
        }
      }
      #if(verbose) cat(tu, "at marker", m , wilcox.test(testQTL[tuID, m], testQTL[exonID[-which(rawexp[exonID, "tu"] == tu)], m])$p.value, "\n")
      resMatrix <- rbind(resMatrix, res1vr, res2vr)#, res1v2
      resRN <- c(resRN, paste0(filename, "_", tu, "_1vr"), paste0(filename, "_", tu, "_2vr"))#, paste0(filename, "_", tu, "_1v2")
    }
  }
  rownames(resMatrix) <- resRN
  return(resMatrix)
}


#QTL
#get the list of QTLGene
toTest = "QTL"; threTest = 8
for(chr in 1:5){
  genenames <- gsub(".txt", "", expGeneList[[chr]])
  QTLGene <- NULL
  resMatrix <- NULL
  for(filename in genenames){
    if(expGeneTest(filename, M = 4, P = 3, N = 2) && potentialQTLGene(filename, N = 2, toTest = "QTL", threTest = 8, Q = 1, S = 3)){
    #TRUE -> it's a portentialQTL/Int gene
    cat(filename, "\n")
    QTLGene <- c(QTLGene, filename)
    resMatrix <- rbind(resMatrix, ttestQTLGene(filename, toTest = "QTL", P = 3))
    }
  }
  #change name when change parameters!!!
  save(QTLGene, file=paste0("Data/countQTL/main", toTest, "_chr", chr, "_fromALL_thre", threTest, ".Rdata"))
  write.table(resMatrix, file = paste0("Data/countQTL/main", toTest, "_chr", chr, "_fromALL_ttest.txt"), sep = " ", row.names = TRUE, col.names = FALSE)
}

toTest = "QTL"; threTest = 8
for(chr in 1:5){
  resMatrix <- NULL
  load(file = paste0("Data/countQTL/main", toTest, "_chr", chr, "_thre", threTest, ".Rdata"))
  
  for(filename in QTLGene){
    resMatrix <- rbind(resMatrix, wilcoxtestQTLGene(filename, toTest = "QTL", P = 3))
    cat(filename, "finished!\n")
  }
  #resMatrix
  write.table(resMatrix, file = paste0("Data/countQTL/main", toTest, "_chr", chr, "_wilcoxExp.txt"), sep = " ", row.names = TRUE, col.names = FALSE)
}


#Int
#get the list of QTLGene
toTest = "Int"; threTest = 11.6
for(chr in 1:5){
  genenames <- gsub(".txt", "", expGeneList[[chr]])
  QTLGene <- NULL
  resMatrix <- NULL
  for(filename in genenames){
    if(expGeneTest(filename, M = 4, P = 3, N = 2) && potentialQTLGene(filename, N = 2, toTest = "Int", threTest = 11.6, Q = 1, S = 3)){
    #TRUE -> it's a portentialQTL/Int gene
    cat(filename, "\n")
    QTLGene <- c(QTLGene, filename)
    resMatrix <- rbind(resMatrix, ttestQTLGene(filename, toTest = "Int", P = 3))
    }
  }
  #change name when change parameters!!!
  save(QTLGene, file=paste0("Data/countQTL/main", toTest, "_chr", chr, "_fromALL_thre", threTest, ".Rdata"))
  write.table(resMatrix, file = paste0("Data/countQTL/main", toTest, "_chr", chr, "_fromALL_ttest.txt"), sep = " ", row.names = TRUE, col.names = FALSE)
}

toTest = "Int"; threTest = 11.6
for(chr in 1:5){
  resMatrix <- NULL
  load(file = paste0("Data/countQTL/main", toTest, "_chr", chr, "_thre", threTest, ".Rdata"))
  
  for(filename in QTLGene){
    resMatrix <- rbind(resMatrix, wilcoxtestQTLGene(filename, toTest = "Int", P = 3))
    cat(filename, "finished!\n")
  }
  #resMatrix
  write.table(resMatrix, file = paste0("Data/countQTL/main", toTest, "_chr", chr, "_wilcoxExp.txt"), sep = " ", row.names = TRUE, col.names = FALSE)
}
