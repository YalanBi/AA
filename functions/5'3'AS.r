#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 19-06-2013
# first written: 21-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#********************************************** This is used to find alternative splicing at 5' or 3' **********************************************#


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

#to find max difference between each probe in exp, for grouping
findMAXdff <- function(toGroup = ind, exp_data = rawexp[ ,17:164], verbose = FALSE){
  dff <- NULL
  for(n in 2:length(toGroup)){
    dff <- c(dff, abs(sum(exp_data[toGroup[n], ] - exp_data[toGroup[n-1], ])))
    if(verbose) cat("difference between p", toGroup[n], "and p", toGroup[n-1], "is", sum(exp_data[toGroup[n], ]) - sum(exp_data[toGroup[n-1], ]), "\n")
  }
  if(verbose) cat("so dff is:", dff, "and max dff is between p", toGroup[which.max(dff)], "and p", toGroup[which.max(dff)+1], "-", max(dff), "\n")
  
  return(which.max(dff))
}
#findMAXdff(toGroup = ind, exp_data = rawexp[ ,17:164], verbose = T)

#5'/3' AS test
#annotation: useForTest = unlist -> use all individuals to do test, better than mean/median
#            whichTest <- wilcox.test/ t.test; one-side!!!
testAS <- function(exp_data = rawexp[ ,17:164], testProbes, restProbes, ind = ind_env, useForTest = unlist, whichTest = wilcox.test, alternative = "less", verbose = FALSE){
  testPart <- apply(exp_data[testProbes, ind], 2, useForTest)
  if(verbose) cat("We are testProbes:", testProbes, "\n")
  restPart <- apply(exp_data[restProbes, ind], 2, useForTest)
  if(verbose) cat("We are restProbes:", restProbes, "\n")
  return(-log10(whichTest(testPart, restPart, alternative) $ p.value))
}
#testAS(exp_data = rawexp[ ,17:164], testProbes, restProbes, ind = ind_env, useForTest = unlist, whichTest = wilcox.test, alternative = "less", verbose = TRUE)

#simple: slimlar with expGeneTestSimple
expExonTestSimple <- function(ind, exp_data = rawexp[ ,17:164], use = median, expThre = 5, verbose = FALSE){
  if(length(ind) > 0){
    expDegree <- use(apply(exp_data[ind, ], 1, unlist))
    if(expDegree >= expThre){
      if(verbose) cat("this exon has p", ind, "and exp degree is", expDegree, ", higher than", expThre, "!\n")
      return(TRUE)
    } else return(FALSE)
  } else return(FALSE)
}
#expExonTestSimple(ind, exp_data = rawexp[ ,17:164], use = median, expThre = 5, verbose = TRUE)

splicingTest <- function(filename, goal, verbose = FALSE, ...){
  chr <- as.numeric(gsub("AT", "", strsplit(filename, "G")[[1]][1]))
  rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  probes_dir <- probesDir(rawexp)
  #if(verbose) cat("We have rightDir probes:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #if(verbose) cat("We have exon probes:", exonID, "\n")
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #if(verbose) cat("We have exons:", uniqueExon, "\n")
  
  #for 5'/3' AS, at least 2 exons in a gene!!!
  if(length(uniqueExon) >= 2){
    if(verbose) cat(filename, "have", length(uniqueExon), "exons, >= 2!\n")
    
    if(goal == "skippingExon"){
      P <- 3
      exonTestRange <- uniqueExon
    }
    if(goal == "35AS"){
      P <- 6
      exonTestRange <- uniqueExon[c(1, length(uniqueExon))]
    }
    if(verbose) cat(goal, "test, among:", exonTestRange, "\n")
    
    resmatrix <- NULL
    rownameList <- NULL
    
    for(testExon in exonTestRange){
      ind <- exonID[rawexp[exonID, "tu"] == testExon]
      #if(verbose) cat(testExon, "has probes", ind, "\n")
      
      #at least 6 probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
      if(length(ind) >= P && expExonTestSimple(ind, exp_data = rawexp[ ,17:164], use = median, expThre = 5, verbose)){
        if(verbose) cat("\t***I'm", testExon, ", has", length(ind), "good probes, separate me and check if i'm possible for later test!\n")
        
        if(goal == "skippingExon"){
          #to find max difference between each probe in exp, to group them
          testProbes <- ind
          if(verbose) cat("I'm testProbes, I have p", testProbes, "\n")
          restProbes <- exonID[!exonID %in% ind]
          if(verbose) cat("I'm restProbes, I have p", restProbes, "\n")
        }
        if(goal == "35AS"){
          separatePoint <- findMAXdff(toGroup = ind, exp_data = rawexp[ ,17:164], verbose)
          testProbes <- ind[1:separatePoint]
          #if(verbose) cat("I'm testpart, I have p", testProbes, "\n")
          restProbes <- ind[-(1:separatePoint)]
          #if(verbose) cat("I'm restpart, I have p", restProbes, "\n")
        }
        
        #>= 3 probes left in each group, remember the gene name and do t.test
        if(length(testProbes) >= 3 && length(restProbes) >= 3){
          separateProbe <- max(testProbes)
          if(verbose) cat(" =>separate at p", separateProbe, ", have >= 3 good probes in each group, ready for test!\n")
          rownameList <- c(rownameList, filename)
          res <- separateProbe
          for(env in 1:4){
            ind_env <- which(as.numeric(menvironment) == env)
            res <- c(res, testAS(exp_data = rawexp[ ,17:164], testProbes, restProbes, ind = ind_env, verbose, ...))
            #if(verbose) cat("env", env, ":", testAS(exp_data = rawexp[ ,17:164], testProbes, restProbes, ind = ind_env, ...), "\n")
          }
          resmatrix <- rbind(resmatrix, res)
        } else if(verbose) cat(" =>after grouping don't have enough probes in each part T^T\n")
      } else if(verbose) cat("\t***I'm", testExon, ", not enough probes/not expressed T^T\n")
    }
    rownames(resmatrix) <- rownameList
    return(resmatrix)
  } else if(verbose) cat("we don't have enough exons T^T\n")
}
#splicingTest(filename, goal = "35AS", useForTest = unlist, whichTest = wilcox.test, alternative = "less", verbose = TRUE)


#*************************************************************** test part ***************************************************************#
#Before start, mean, median or all individuals, change!!!
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  #filename = "AT1G01010"
  for(filename in genenames){
    res <- splicingTest(filename, goal = "35AS", useForTest = unlist, whichTest = wilcox.test, alternative = "less")
    if(!is.null(res)){
      resmatrix <- rbind(resmatrix, res)
      cat(filename, "is tested\n")
    }
  }
  colnames(resmatrix) <- c("sepProbe", "6H", "Dry_AR", "Dry_Fresh", "RP")
  
  write.table(resmatrix, file=paste0("Data/53terminalAS/53AS_chr", chr, "_wt_less.txt"), sep="\t") #********** change!!! **********#
  
  et <- proc.time()[3]
  cat("chr", chr, "finished in", et-st, "s\n\n")
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
