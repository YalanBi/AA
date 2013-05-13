#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 13-05-2013
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

#cassette exon test
testCassette <- function(expdata = rawexp[exonID, ind_env + 16], ind, use){
  testProbe <- apply(expdata[ind, ], 2, use)
  #cat(testProbe, "\n")
  otherProbe <- apply(expdata[!ind, ], 2, use)
  #cat(otherProbe, "\n")
  return(-log10(t.test(testProbe, otherProbe, alternative="less")$p.value))
}


#*************************************************************** load part ***************************************************************#
#load exp genes
load(file="Data/fullModeMapping/expGenes.Rdata")


#*************************************************************** test part ***************************************************************#
#number of t.test we have done, used for FDR
cnt <- c(0, 0, 0, 0, 0)

#Before start, mean, median or all individuals, change!!!
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrix <- NULL
  rownameList <- NULL
  
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
    
    #at least 3 exons in a gene!!!
    if(length(uniqueExon) >= 3){
      #cat("we have >= 3 exons!\n")
      
      for(exon in uniqueExon[2: (length(uniqueExon)-1)]){
        #ind <- judge which probe in exonID is of current exon name (T/F)
        ind <- rawexp[exonID, "tu"] == exon
        #cat(as.character(exon), "has probes", exonID[ind], "\n")
        
        #at least 3 probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
        if(length(which(ind)) >= 3){
          #cat("I'm", exon, ">= 3 good probes, t.test me and remember my tu name!\n")
          rownameList <- c(rownameList, paste0(gsub(".txt", "", filename), "_", exon))
          
          res <- NULL
          for(env in 1:4){
            ind_env <- which(as.numeric(menvironment) == env)
            
            #use mean/median/all individuals(unlist) to do the t.test for cassette exon
            res <- c(res, testCassette(expdata = rawexp[exonID, ind_env + 16], ind, use = unlist)) #********** change!!! **********#
            
            #t.test once, counter plus 1!!! to calculate the number of t.test we have done. for FDR
            cnt[chr] <- cnt[chr] + 1
          }
          
        resmatrix <- rbind(resmatrix, res)
        }
      }
    }
  }
  
  rownames(resmatrix) <- rownameList
  colnames(resmatrix) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
  write.table(resmatrix, file=paste0("Data/cassetteExon/cassetteExon_chr", chr, "_allind.txt"), sep="\t") #********** change!!! **********#
  
  et <- proc.time()[3]
  cat("and finished in", et-st, "s\n")
}

cnt #cnt = length(rownameList * 4)
[1] 9344 5260 7428 5568 8416
#sum(cnt)=36016


#count nTestTus and nTestGenes
ce_threshold = 5.86

nTestTus <- NULL
nTestGenes <- NULL

matrixSigGenes <- NULL
matrixSigTus <- NULL

for(chr in 1:5){
  cematrix <- read.table(paste0("Data/cassetteExon/cassetteExon_chr", chr, "_allind.txt"), row.names=1, header=T) #********** change!!! **********#
  
  nTestTus <- c(nTestTus, nrow(cematrix))
  nTestGenes <- c(nTestGenes, length(unique(unlist(lapply(strsplit(rownames(cematrix),"_"),"[[",1)))))
  
  #for all genes which have cassette exons in one env
  nSigGenes <- NULL
  nSigTus <- NULL
  #for all genes which have cassette exons in one env
  for(env in 1:4){
    nSigGenes <- c(nSigGenes, length(unique(unlist(lapply(strsplit(rownames(cematrix)[which(cematrix[ ,env] >= ce_threshold)],"_"),"[[",1)))))
    nSigTus <- c(nSigTus, length(which(cematrix[ ,env] >= ce_threshold)))
  }
  matrixSigGenes <- rbind(matrixSigGenes, nSigGenes)
  matrixSigTus <- rbind(matrixSigTus, nSigTus)
}

nTestTus
[1] 2336 1315 1857 1392 2104

nTestGenes
[1] 1338  742 1058  784 1227

matrixSigGenes(AS genes)
      Env1  Env2  Env3  Env4
chr1  488   526   530   528
chr2  302   318   322   296
chr3  386   413   433   395
chr4  291   304   309   302
chr5  447   483   475   453

matrixSigTus(AS exons)
      Env1  Env2  Env3  Env4
chr1  606   647   647   657
chr2  356   378   381   357
chr3  466   501   524   469
chr4  355   380   383   358
chr5  540   586   569   536



#count ngene having cassette exon
geneCExon_cnt <- NULL
for(chr in 1:5){
  cematrix <- read.table(paste0("Data/cassetteExon/cassetteExon_chr", chr, "_allind.txt"), row.names=1, header=T) #********** change!!! **********#
  
  
  geneCExon_cnt <- rbind(geneCExon_cnt, nSigGenes)
  
}
geneCExon_cnt








#************************************************************* plot exp part *************************************************************#
plotcExonExp <- function(chr, filename, ce_threshold, use = "median"){
  rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  newexp <- rawexp[,17:164]
  
  probes_dir <- probesDir(rawexp)
  #cat(filename, "\nprobeDir:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
  #cat("exons of right direction:", exonID, "\n")
  
  #uniqueExon <- all tu names of exon probes
  uniqueExon <- unique(rawexp[exonID, "tu"])
  #cat("tu names:", as.character(uniqueExon), "\n")
  
  ind_tu <- grep("tu", rawexp[,"tu"])
  nprobes <- nrow(rawexp)
  
  png(filename = paste0("Data/cassetteExon/chr", chr, "/", filename,"_cExonExp_", use, ce_threshold, ".png"), width = 1024, height = 1024, bg = "white")
  plot(c(0.5, nprobes+0.5), c(min(newexp)-0.2, max(newexp))+0.2, xaxt='n', xlab="Probes", ylab="Intensity", cex.axis=1, cex.lab=1.5, cex.main=2, las=1, mgp=c(3,1,0), tck=-0.017, t="n", main=filename)
  
  for(p in 1:nprobes){
    #background for introns
    if(!p %in% ind_tu){
      rect((p-0.5), -3, (p+0.5), max(newexp)* 2, col=grey(0.85), border = "transparent")
    }
    if(p %in% probes_dir){
      points(rep(p,148)+0.06*as.numeric(menvironment)-05, newexp[p,], t='p', col=menvironment, pch=20, cex=0.7)
    }
  }
  for(exon in uniqueExon){
    #ind <- judge which probe in exonID is of current exon name (T/F)
    ind <- rawexp[exonID, "tu"] == exon
    #cat(as.character(exon), "has probes", exonID[ind], "\n")
    
    for(env in 1:4){
      ind_env <- which(as.numeric(menvironment) == env)
      
      #use mean/median to test cassette
      lines(c(min(which(rawexp[,"tu"] == exon))-0.5+0*env, max(which(rawexp[,"tu"] == exon))+0*env), c(mean(unlist(newexp[exonID[ind], ind_env])), mean(unlist(newexp[exonID[ind], ind_env]))), col=env, lwd=2)       
      
      cExoninEnv <- unlist(lapply(strsplit(row.names(cematrix[which(cematrix[,env]>=ce_threshold),])[grepl(filename, row.names(cematrix[which(cematrix[,env]>=ce_threshold),]))],"_"),"[[",2))
      if(exon %in% cExoninEnv){
        text(x=median(exonID[ind]), y=max(newexp)-0*env, labels=exon, col=env, cex=1)
      }
    }
  }
  axis(1, at=probes_dir, labels=row.names(rawexp)[probes_dir], cex.axis=1, las=2, tck=0.005)
  box()
  dev.off()
}


use = "median"
for(chr in 1:5){
  cematrix <- read.table(paste0("Data/cassetteExon/cassetteExon_chr", chr, "_", use, "_D2p.txt"), row.names=1, header=T)
  
  #for all genes which have cassette exons in at least one env
  plotGenenames <- unique(unlist(lapply(strsplit(rownames(cematrix)[which(cematrix >= ce_threshold, arr.ind=T)],"_"),"[[",1)))
  #for genes having cassette exons and the -log10(P) is round the ce_threshold
  #ce_threshold = 6.3
  #plotGenenames <- unique(unlist(lapply(strsplit(rownames(which(round(cematrix, digits=1)==ce_threshold, arr.ind=T)),"_"),"[[",1)))
  #for genes having cassette exons and the -log10(P) is higher than 55
  #plotGenenames <- unique(unlist(lapply(strsplit(rownames(which(cematrix > 55, arr.ind=T)),"_"),"[[",1)))
  
  #filename(AT1G01010)
  for(filename in plotGenenames[1:15]){
    plotcExonExp(chr, filename, ce_threshold = 6.3, use = "median")
  }
}
