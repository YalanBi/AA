#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 29-05-2013
# first written: 21-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[ ,2]


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


#5'/3' AS test
#t.test; use = unlist -> use all individuals to do t.test, better than mean/median
testASinExon <- function(expdata = rawexp[ ,ind_env + 16], part1, part2, use){
  part1Probes <- apply(expdata[part1, ], 2, use)
  #cat("We are part1Probes:", part1, "\n")
  part2Probes <- apply(expdata[part2, ], 2, use)
  #cat("We are part2Probes:", part2, "\n")
  return(-log10(t.test(part1Probes, part2Probes) $ p.value))
}
#wilcox.test; use = unlist -> use all individuals to do t.test, better than mean/median
testASinExon <- function(expdata = rawexp[ ,ind_env + 16], part1, part2, use){
  part1Probes <- apply(expdata[part1, ], 2, use)
  #cat("We are part1Probes:", part1, "\n")
  part2Probes <- apply(expdata[part2, ], 2, use)
  #cat("We are part2Probes:", part2, "\n")
  return(-log10(wilcox.test(part1Probes, part2Probes) $ p.value))
}


#*************************************************************** load part ***************************************************************#
#load exp genes
load(file="Data/fullModeMapping/expGenes.Rdata")


#*************************************************************** test part ***************************************************************#
#Before start, mean, median or all individuals, change!!!
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  rownameList <- NULL
  
  #filename = "AT1G01010.txt"
  for(filename in genenames){
    #cat(filename, "...\n")
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename), row.names=1, header=T)
    #cat("rawexp loading succeed!\n")
    
    #get probes' and exons' ID
    probes_dir <- probesDir(rawexp)
    #cat("probes of right direction:", probes_dir, "\n")
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
    #cat("exons of right direction:", exonID, "\n")
    
    #uniqueExon <- all tu names of exon probes
    uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
    #cat(" =>tu names:", as.character(uniqueExon), "\n")
    
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
