#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 18-06-2013
# first written: 18-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#********************************************** This is used to find genetic regulated AS at 5' or 3' **********************************************#


setwd("D:/Arabidopsis Arrays")

menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]

#load exp genes
load(file="Data/ExpGenes/expGenes_simple.Rdata")

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

findMAXdff <- function(probesToGroup = ind, exp_data = rawexp[ ,17:164], verbose = FALSE){
  #to find max difference between each probe in exp, to group them
  dff <- NULL
  for(n in 2:length(probesToGroup)){
    dff <- c(dff, abs(sum(exp_data[probesToGroup[n], ] - exp_data[probesToGroup[n-1], ])))
    if(verbose) cat("difference between p", probesToGroup[n], "and p", probesToGroup[n-1], "is", sum(exp_data[probesToGroup[n], ]) - sum(exp_data[probesToGroup[n-1], ]), "\n")
  }
  if(verbose) cat("so dff is:", dff, "and max dff is between p", probesToGroup[which.max(dff)], "and p", probesToGroup[which.max(dff)+1], "-", max(dff), "\n")
  
  return(which.max(dff))
}
#findMAXdff(probesToGroup = ind, exp_data = rawexp[ ,17:164], verbose = T)

#5'/3' AS test
#annotation: useForTest = unlist -> use all individuals to do test, better than mean/median
#            whichTest <- wilcox.test/ t.test; one-side!!!
testAS <- function(rawexp, testProbes, restProbes, ind = ind_env, useForTest = unlist, whichTest = wilcox.test, verbose = FALSE){
  testPart <- apply(rawexp[testProbes, ind+16], 2, useForTest)
  if(verbose) cat("We are testProbes:", testProbes, "\n")
  restPart <- apply(rawexp[restProbes, ind+16], 2, useForTest)
  if(verbose) cat("We are restProbes:", restProbes, "\n")
  return(-log10(whichTest(testPart, restPart, alternative = "less") $ p.value))
}


function(filename, N, exonTestRange, ){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  probes_dir <- probesDir(rawexp)
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  
  if(goal == "skippingExon"){
    P <- 3
    exonTestRange <- uniqueExon
  }
  if(goal == "35AS"){
    P <- 6
    exonTestRange <- uniqueExon[c(1, length(uniqueExon))]
  }
  
  #for 5'/3' AS, at least 2 exons in a gene!!!
  if(length(uniqueExon) >= 2){
    #cat(" we have", length(uniqueExon), "exons, >= 2!\n")
    
    for(testExon in exonTestRange){
      ind <- exonID[rawexp[exonID, "tu"] == testExon]
      if(verbose) cat(testExon, "has probes", ind, "\n")
      
      #at least 6 probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
      if(length(ind) >= P){
        if(verbose) cat("\t***I'm", testExon, ", has", length(ind), "good probes, separate me and check if i'm possible for later test!\n")
        
        #to find max difference between each probe in exp, to group them
        separateProbe <- findMAXdff(probesToGroup = ind, exp_data = rawexp[ ,17:164])
        part1Probes <- ind[1:separateProbe]
        if(verbose) cat("I'm part1, I have", part1Probes, "\n")
        part2Probes <- ind[-(1:separateProbe)]
        if(verbose) cat("I'm part2, I have", part2Probes, "\n")
        
        #>= 3 probes left in each group, remember the gene name and do t.test
        if(length(part1) >= 3 && length(part2) >= 3){
          cat(" =>I'm", filename, testExon, ", I have >= 3 good probes in each group, t.test me and remember my gene_tu name!\n")
          rownameList <- c(rownameList, filename)
          
          res <- separateProbe
          for(env in 1:4){
            ind_env <- which(as.numeric(menvironment) == env)
            #cat("Now is env", env, "\n")
            res <- c(res, testAS(rawexp, testProbes, restProbes, ind = ind_env))
          }
          resmatrix <- rbind(resmatrix, res)
        } #else cat(" =>I'm", testExon, ", but after grouping don't have enough probes in each part T^T\n")
      } #else cat("\t***I'm", as.character(testExon), "but not enough probes T^T\n")
    }
  } #else cat("we don't have enough exons T^T\n")
}



#*************************************************************** test part ***************************************************************#
#Before start, mean, median or all individuals, change!!!
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  rownameList <- NULL
  
  #filename = "AT1G01010"
  for(filename in genenames){
    #cat(filename, "starts...\n")
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
    
    probes_dir <- probesDir(rawexp)
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
    uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
    
    #for 5'/3' AS, at least 2 exons in a gene!!!
    if(length(uniqueExon) >= 2){
      #cat(" we have", length(uniqueExon), "exons, >= 2!\n")
      
      for(testExon in uniqueExon[c(1, length(uniqueExon))]){
        #ind <- judge which probe in exonID is of 1st exon name (T/F)
        ind <- rawexp[exonID, "tu"] == testExon
        #cat("", as.character(testExon), "has probes", exonID[ind], "\n")
        
        #at least 6 probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
        if(length(which(ind)) >= 6){
          #cat("\t***I'm", as.character(testExon), ", has >= 6 good probes, separate me to check whether t.test or not!\n")
          
          #to find max difference between each probe in exp, to group them
          dff <- NULL
          for(n in 2:length(exonID[ind])){
          
            #*********************************************** should it be sum(abs(differences)) or abs(sum(differences))? ? ? ***********************************************#
            dff <- c(dff, sum(abs(rawexp[exonID[ind][n], 17:164] - rawexp[exonID[ind][n-1], 17:164])))
            #cat("  difference between p", exonID[ind][n], "and p", exonID[ind][n-1], "is", sum(rawexp[exonID[ind][n], 17:164]) - sum(rawexp[exonID[ind][n-1], 17:164]), "\n")
          }
          #cat(" so dff is:", dff, "\n")
          
          #find the max difference and separate into 2 parts
          #cat("  max dff is between p", exonID[ind][which.max(dff)], "and p", exonID[ind][which.max(dff)+1], ", is", max(dff), "\n")
          part1 <- exonID[ind][1:which.max(dff)]
          #cat("  I'm part1, I have", part1, "\n")
          part2 <- exonID[ind][-(1:which.max(dff))]
          #cat("  I'm part2, I have", part2, "\n")
          
          #>= 3 probes left in each group, remember the gene name and do t.test
          if(length(part1) >= 3 && length(part2) >= 3){
            cat(" =>I'm", filename, testExon, ", I have >= 3 good probes in each group, t.test me and remember my gene_tu name!\n")
            rownameList <- c(rownameList, paste0(gsub(".txt", "", filename), "_", match(testExon, uniqueExon)))
            
            res <- NULL
            for(env in 1:4){
              ind_env <- which(as.numeric(menvironment) == env)
              #cat("Now is env", env, "\n")
              res <- c(res, testASinExon(expdata = rawexp[ ,ind_env + 16], part1, part2, use = unlist))
            }
            resmatrix <- rbind(resmatrix, res)
          } #else cat(" =>I'm", testExon, ", but after grouping don't have enough probes in each part T^T\n")
        } #else cat("\t***I'm", as.character(testExon), "but not enough probes T^T\n")
      }
    } #else cat("we don't have enough exons T^T\n")
  }
  rownames(resmatrix) <- rownameList
  colnames(resmatrix) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
  #write.table(resmatrix, file=paste0("Data/53terminalAS/53terminalAS_chr", chr, "_ttest.txt"), sep="\t") #********** change!!! **********#
  write.table(resmatrix, file=paste0("Data/53terminalAS/53terminalAS_chr", chr, "_wtest.txt"), sep="\t") #********** change!!! **********#
  
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
}





#select expressed exon probes
chooseExpExon <- function(exonID, newexp, ind_env, cutoff){
  expExonID <- NULL
  for(p in exonID){
    #cat(p, mean(unlist(newexp[p, ind_env])), "\n")
    
    #cutoff <- cutoff used to check whether mean(each exon of right direction) is high enough to be regarded as expressed or not
    #           then remove low expressed exons
    if(mean(unlist(newexp[p, ind_env])) >= cutoff){
      expExonID <- c(expExonID, p)
    }
  }
  return(expExonID)
}







#*************************************************************** count part ***************************************************************#
nTest <- NULL
nTest3or5Genes <- NULL
nTest5Genes <- NULL
nTest3Genes <- NULL
for(chr in 1:5){
  #aa <- read.table(paste0("Data/53terminalAS/53terminalAS_chr", chr, "_ttest.txt"), row.names=1, header=T)
  aa <- read.table(paste0("Data/53terminalAS/53terminalAS_chr", chr, "_wtest.txt"), row.names=1, header=T)
  nTest <- c(nTest, nrow(aa) * ncol(aa))
  
  rowGenes <- unlist(lapply(strsplit(rownames(aa), "_"), "[[", 1))
  nTest3or5Genes <- c(nTest3or5Genes, length(unique(rowGenes)))
  
  firstGenes <- rowGenes[unlist(lapply(strsplit(rownames(aa), "_"), "[[", 2)) == 1]
  nTest5Genes <- c(nTest5Genes, length(firstGenes))
  
  lastGenes <- rowGenes[unlist(lapply(strsplit(rownames(aa), "_"), "[[", 2)) != 1]
  nTest3Genes <- c(nTest3Genes, length(lastGenes))
}
nTest
[1] 1476  812 1084  712 1252
sum(nTest) = 5336

ps_threshold = -log10(0.05/sum(nTest)) = 5.03

nTest3or5Genes
[1] 355 196 265 174 297
nTest5Genes
[1] 182 100 135  94 160
nTest3Genes
[1] 187 103 136  84 153


matrixSig3or5 <- NULL
matrixSig5 <- NULL
matrixSig3 <- NULL
for(chr in 1:5){
  #aa <- read.table(paste0("Data/53terminalAS/53terminalAS_chr", chr, "_ttest.txt"), row.names=1, header=T)
  aa <- read.table(paste0("Data/53terminalAS/53terminalAS_chr", chr, "_wtest.txt"), row.names=1, header=T)
  
  #for all genes which have cassette exons in one env
  nSig3or5Genes <- NULL
  nSig5Genes <- NULL
  nSig3Genes <- NULL
  #for all genes which have cassette exons in one env
  for(env in 1:4){
    rowSigGenes <- unlist(lapply(strsplit(rownames(aa)[aa[ ,env] >= ps_threshold], "_"), "[[", 1))
    nSig3or5Genes <- c(nSig3or5Genes, length(unique(rowSigGenes)))
    
    firstSigGenes <- rowSigGenes[unlist(lapply(strsplit(rownames(aa)[aa[ ,env] >= ps_threshold], "_"), "[[", 2)) == 1]
    nSig5Genes <- c(nSig5Genes, length(firstSigGenes))
    
    lastSigGenes <- rowSigGenes[unlist(lapply(strsplit(rownames(aa)[aa[ ,env] >= ps_threshold], "_"), "[[", 2)) != 1]
    nSig3Genes <- c(nSig3Genes, length(lastSigGenes))
    
  }
  matrixSig3or5 <- rbind(matrixSig3or5, nSig3or5Genes)
  matrixSig5 <- rbind(matrixSig5, nSig5Genes)
  matrixSig3 <- rbind(matrixSig3, nSig3Genes)
}
rownames(matrixSig3or5) <- c("chr1", "chr2", "chr3", "chr4", "chr5")
colnames(matrixSig3or5) <- c("Env1", "Env2", "Env3", "Env4")
rownames(matrixSig5) <- c("chr1", "chr2", "chr3", "chr4", "chr5")
colnames(matrixSig5) <- c("Env1", "Env2", "Env3", "Env4")
rownames(matrixSig3) <- c("chr1", "chr2", "chr3", "chr4", "chr5")
colnames(matrixSig3) <- c("Env1", "Env2", "Env3", "Env4")

#results of t.test
matrixSig3or5(ngenes having AS at 3|5 site)
      Env1 Env2 Env3 Env4
chr1  242  247  250  234
chr2  131  138  141  133
chr3  170  177  175  166
chr4  120  118  118  110
chr5  208  212  223  196
matrixSig5(ngenes having AS at 5 site)
      Env1 Env2 Env3 Env4
chr1  132  129  131  130
chr2   63   67   67   60
chr3   94   89   89   89
chr4   66   67   66   61
chr5  113  111  118  105
matrixSig3(ngenes having AS at 3 site)
      Env1 Env2 Env3 Env4
chr1  115  124  125  109
chr2   70   73   76   73
chr3   77   89   87   78
chr4   55   53   54   50
chr5  100  106  110   98






#*************************************************************** main idea ***************************************************************#
rr <- NULL
for(x in 1:(nrow(rawexp)-1)){
  rr <- c(rr, sum(abs(rawexp[x,17:164]-rawexp[(x+1),17:164])))
}

#minimum of 4 probes
#Test how to split (using highest difference)
#Split into groups
#test if every group has 2 probes
#YES -> T-Test against the first group (5" start) all the other probes in the gene
#NO -> continue
