#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 16-05-2013
# first written: 02-05-2013
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

#5'/3' exon skipping test
#use = unlist -> use all individuals to do t.test, better than mean/median
testSkipping <- function(expdata = rawexp[exonID,ind_env + 16], ind, use){
  testProbes <- apply(expdata[ind, ], 2, use)
  #cat("We are testProbes:", exonID[ind], "\n")
  otherProbes <- apply(expdata[!ind, ], 2, use)
  #cat("We are otherProbes:", exonID[!ind], "\n")
  return(-log10(t.test(testProbes, otherProbes, alternative="less") $ p.value))
}


#*************************************************************** load part ***************************************************************#
#load exp genes
load(file="Data/fullModeMapping/expGenes.Rdata")


#*************************************************************** test part ***************************************************************#
#number of t.test we have done, used for FDR
cntF <- c(0, 0, 0, 0, 0)
cntL <- c(0, 0, 0, 0, 0)

#Before start, mean, median or all individuals, change!!!
for(chr in 1:5){
  st <- proc.time()[3]
  cat("chr", chr, "starts...\n")
  
  genenames <- expGeneList[[chr]]
  resmatrixF <- NULL
  rownameListF <- NULL
  resmatrixL <- NULL
  rownameListL <- NULL
  
  for(filename in genenames){
    cat(filename, "...\n")
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename), row.names=1, header=T)
    #cat("rawexp loading succeed!\n")
    
    #get probes' and exons' ID
    probes_dir <- probesDir(rawexp)
    #cat("probes of right direction:", probes_dir, "\n")
    exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
    #cat("exons of right direction:", exonID, "\n")
    
    #uniqueExon <- all tu names of exon probes
    uniqueExon <- unique(rawexp[exonID,"tu"])
    #cat("tu names:", as.character(uniqueExon), "\n")
    
    #for skipping exon, at least 2 exons in a gene!!!
    if(length(uniqueExon) >= 2){
      #cat("we have >= 2 exons!\n")
      
    #*************************first, test the 1st exon/5' exon*************************#
      #indF <- judge which probe in exonID is of 1st exon name (T/F)
      indF <- rawexp[exonID, "tu"] == uniqueExon[1]
      #cat(as.character(uniqueExon[1]), "has probes", exonID[indF], "\n")
      
      #at least 3 probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
      if(length(which(indF)) >= 3){
        #cat("I'm", as.character(uniqueExon[1]), ">= 3 good probes, t.test me and remember my gene name!\n")
        rownameListF <- c(rownameListF, gsub(".txt", "", filename))
        
        resF <- NULL
        for(env in 1:4){
          ind_env <- which(as.numeric(menvironment) == env)
          
          #use mean/median/all individuals(unlist) to do the t.test for cassette exon
          resF <- c(resF, testSkipping(expdata = rawexp[exonID,(ind_env + 16)], ind = indF, use = unlist)) #********** change!!! **********#
          
          #t.test once, counter plus 1!!! to calculate the number of t.test we have done. for FDR
          cntF[chr] <- cntF[chr] + 1
        }
      resmatrixF <- rbind(resmatrixF, resF)
      }
      
    #*************************next, test the last exon/3' exon*************************#
      #indL <- judge which probe in exonID is of last exon name (T/F)
      indL <- rawexp[exonID, "tu"] == uniqueExon[length(uniqueExon)]
      #cat(as.character(uniqueExon[length(uniqueExon)]), "has probes", exonID[indL], "\n")
      
      #at least 3 probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
      if(length(which(indL)) >= 3){
        #cat("I'm", as.character(uniqueExon[length(uniqueExon)]), ">= 3 good probes, t.test me and remember my gene name!\n")
        rownameListL <- c(rownameListL, gsub(".txt", "", filename))
        
        resL <- NULL
        for(env in 1:4){
          ind_env <- which(as.numeric(menvironment) == env)
          
          #use mean/median/all individuals(unlist) to do the t.test for cassette exon
          resL <- c(resL, testSkipping(expdata = rawexp[exonID,(ind_env + 16)], ind = indL, use = unlist)) #********** change!!! **********#
          
          #t.test once, counter plus 1!!! to calculate the number of t.test we have done. for FDR
          cntL[chr] <- cntL[chr] + 1
        }
      resmatrixL <- rbind(resmatrixL, resL)
      }
    }
  }
  
  rownames(resmatrixF) <- rownameListF
  colnames(resmatrixF) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
  write.table(resmatrixF, file=paste0("Data/skippingExon/5skippingExon_chr", chr, "_allind.txt"), sep="\t") #********** change!!! **********#
  
  rownames(resmatrixL) <- rownameListL
  colnames(resmatrixL) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
  write.table(resmatrixL, file=paste0("Data/skippingExon/3skippingExon_chr", chr, "_allind.txt"), sep="\t") #********** change!!! **********#
  
  et <- proc.time()[3]
  cat("and finished in", et-st, "s\n")
}

cntF #cnt = length(rownameList) * 4
[1] 6252 3740 4876 3640 5696
#sum(cntF)=24204
cntL
[1] 6484 3652 5028 3968 5936
#sum(cntF)=25068


#************************************************************* analysis part *************************************************************#
#count nSkipF and nSkipL
fs_threshold = 5.68
ls_threshold = 5.70

nTestF <- NULL
nTestL <- NULL

matrixSkipF <- NULL
matrixSkipL <- NULL

for(chr in 1:5){
  fsmatrix <- read.table(paste0("Data/skippingExon/5skippingExon_chr", chr, "_allind.txt"), row.names=1, header=T) #********** change!!! **********#
  lsmatrix <- read.table(paste0("Data/skippingExon/3skippingExon_chr", chr, "_allind.txt"), row.names=1, header=T) #********** change!!! **********#
  
  nTestF <- c(nTestF, nrow(fsmatrix))
  nTestL <- c(nTestL, nrow(lsmatrix))
  
  #for all genes which have cassette exons in one env
  nSkipF <- NULL
  nSkipL <- NULL
  #for all genes which have cassette exons in one env
  for(env in 1:4){
    nSkipF <- c(nSkipF, length(which(fsmatrix[ ,env] >= fs_threshold)))
    nSkipL <- c(nSkipL, length(which(lsmatrix[ ,env] >= ls_threshold)))
  }
  matrixSkipF <- rbind(matrixSkipF, nSkipF)
  matrixSkipL <- rbind(matrixSkipL, nSkipL)
}

nTestF
[1] 1563  935 1219  910 1424

nTestL
[1] 1621  913 1257  992 1484

matrixSkipF#(5' skipping genes)
      chr1  chr2  chr3  chr4  chr5
Env1  312   198   260   199   320
Env2  321   199   271   198   326
Env3  321   212   271   206   331
Env4  310   201   251   203   325

matrixSkipL#(3' skipping genes)
      chr1  chr2  chr3  chr4  chr5
Env1  343   196   250   211   304
Env2  335   196   255   228   309
Env3  338   195   259   232   315
Env4  358   190   251   212   309




#********************************************************** merge 5'/3' together **********************************************************#
fs_threshold = 5.68
ls_threshold = 5.70
nMergeTest <- NULL
for(chr in 1:5){
  fsmatrix <- read.table(paste0("Data/skippingExon/5skippingExon_chr", chr, "_allind.txt"), row.names=1, header=T) #********** change!!! **********#
  lsmatrix <- read.table(paste0("Data/skippingExon/3skippingExon_chr", chr, "_allind.txt"), row.names=1, header=T) #********** change!!! **********#
  
  frn <- rownames(fsmatrix)
  lrn <- rownames(lsmatrix)
  namerange <- unique(c(frn, lrn))[order(unique(c(frn, lrn)))]#sort()????????????
  nMergeTest <- c(nMergeTest, length(unique(c(frn, lrn))))
  
  mergeSigmatrix <- NULL
  mergeSigRnames <- NULL
  for(g in namerange){
    if(g %in% frn && any(fsmatrix[g, ] >= fs_threshold)){
      mergeSigmatrix <- rbind(mergeSigmatrix, fsmatrix[g, ])
      mergeSigRnames <- c(mergeSigRnames, paste0(g, "_first"))
    }
    if(g %in% lrn && any(lsmatrix[g, ] >= ls_threshold)){
      mergeSigmatrix <- rbind(mergeSigmatrix, lsmatrix[g, ])
      mergeSigRnames <- c(mergeSigRnames, paste0(g, "_last"))
    }
  }
  rownames(mergeSigmatrix) <- mergeSigRnames
  colnames(mergeSigmatrix) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
  write.table(mergeSigmatrix, file=paste0("Data/skippingExon/skippingExon_chr", chr, "_merge.txt"), sep="\t")
}
nMergeTest
[1] 2270 1307 1784 1351 2056

totmatrix <- NULL
for(chr in 1:5){
  mmatrix <- read.table(paste0("Data/skippingExon/skippingExon_chr", chr, "_merge.txt"), row.names=1, header=T)
  ntot <- NULL
  for(env in 1:4){
    ntot <- c(ntot, length(unique(c(unlist(lapply(strsplit(rownames(mmatrix)[grepl("f", rownames(mmatrix))][mmatrix[grepl("f", rownames(mmatrix)),env] >= fs_threshold], "_"), "[[", 1)), unlist(lapply(strsplit(rownames(mmatrix)[grepl("l", rownames(mmatrix))][mmatrix[grepl("l", rownames(mmatrix)),1] >= ls_threshold], "_"), "[[", 1))))))
  }
  totmatrix <- rbind(totmatrix, ntot)
}
totmatrix
#5'or3'skipping genes
     Env1 Env2 Env3 Env4
chr1  646  651  651  639
chr2  386  388  400  388
chr3  506  516  516  494
chr4  395  394  403  398
chr5  615  621  626  618




#************************************************************* plot exp part *************************************************************#
#plot 4 env in 1 plot
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


plotcExonExp <- function(chr, filename, fs_threshold = 5.68, ls_threshold = 5.70){
  rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  newexp <- rawexp[ ,17:164]
  
  probes_dir <- probesDir(rawexp)
  #cat(filename, "\nprobeDir:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #cat("exons of right direction:", exonID, "\n")
  
  #uniqueExon <- all tu names of exon probes
  uniqueExon <- unique(rawexp[exonID,"tu"])
  #cat("tu names:", as.character(uniqueExon), "\n")
  
  ind_tu <- grep("tu", rawexp[ ,"tu"])
  nprobes <- nrow(rawexp)
  
  png(filename = paste0("Data/skippingExon/chr", chr, "/", filename,"_seExp_allind.png"), width = 1024, height = 1024, bg = "white")
  par(oma=c(0.25, 0.5, 0, 0.25))
  plot(c(0.5, nprobes+0.5), c(min(newexp)-0.2, max(newexp)+0.2), xaxt = 'n', xlab = "Probes", ylab = "Intensity", cex.axis = 1, cex.lab = 1.5, cex.main = 2, las = 1, mgp = c(3,0.75,0), tck = -0.017, t = "n", main = filename)
  for(p in 1:nprobes){
    #background for introns
    if(!p %in% ind_tu){
      rect((p-0.5), -3, (p+0.5), max(newexp) * 2, col = grey(0.85), border = "transparent")
    }
    if(p %in% probes_dir){
      points(rep(p, 148)+0.06*as.numeric(menvironment)-0.15, newexp[p, ], t = 'p', col = menvironment, pch = 20, cex = 0.7)
    }
  }
  
  for(exon in uniqueExon){
    #ind <- judge which probe in exonID is of current exon name (T/F)
    ind <- rawexp[exonID,"tu"] == exon
    #cat(as.character(exon), "has probes", exonID[ind], "\n")
    if(length(exonID[ind]) >= 3 && (exon == uniqueExon[1] || exon == uniqueExon[length(uniqueExon)])) text(x = median(which(rawexp[ ,"tu"] == exon)), y = max(newexp) + 0.15, labels = exon, col = "magenta4", cex = 1)
    for(env in 1:4){
      ind_env <- which(as.numeric(menvironment) == env)
      
      #use mean/median to test cassette
      lines(c(min(which(rawexp[ ,"tu"] == exon)) - 0.58 + 0.08 * env, max(which(rawexp[ ,"tu"] == exon)) + 0.18 + 0.08*env), c(mean(unlist(newexp[exonID[ind],ind_env])), mean(unlist(newexp[exonID[ind],ind_env]))), col = env, lwd = 2)       
      
      if(exon == uniqueExon[1] && filename %in% rownames(fsmatrix) && fsmatrix[filename,env] >= fs_threshold) text(x = median(which(rawexp[ ,"tu"] == exon)), y = max(newexp) - 0.15 * env, labels = round(fsmatrix[filename,env], digits = 1), col = env, cex = 0.8)
      if(exon == uniqueExon[length(uniqueExon)] && filename %in% rownames(lsmatrix) && lsmatrix[filename,env] >= ls_threshold) text(x = median(which(rawexp[ ,"tu"] == exon)), y = max(newexp) - 0.15 * env, labels = round(lsmatrix[filename,env], digits = 1), col = env, cex = 0.8)
    }
  }
  axis(1, at = probes_dir, labels = row.names(rawexp)[probes_dir], cex.axis = 1, las = 2, tck = 0.005)
  box()
  dev.off()
}


fs_threshold = 5.68
ls_threshold = 5.70
for(chr in 1:1){
  fsmatrix <- read.table(paste0("Data/skippingExon/5skippingExon_chr", chr, "_allind.txt"), row.names=1, header=T)
  lsmatrix <- read.table(paste0("Data/skippingExon/3skippingExon_chr", chr, "_allind.txt"), row.names=1, header=T)
  
  #for all genes which have cassette exons in at least one env
  skipGenes <- unique(c(rownames(which(fsmatrix >= fs_threshold, arr.ind=TRUE)), rownames(which(lsmatrix >= ls_threshold, arr.ind=TRUE))))
  plotGenenames <- skipGenes[order(skipGenes)]
  #for genes having cassette exons and the -log10(P) is round the ce_threshold
  #aroundSigGenes <- unique(c(rownames(which(round(fsmatrix, digits=2) >= fs_threshold & round(fsmatrix, digits=2) < fs_threshold + 0.2, arr.ind=TRUE)), rownames(which(round(lsmatrix, digits=2) >= ls_threshold & round(lsmatrix, digits=2) < ls_threshold + 0.2, arr.ind=TRUE))))
  #plotGenenames <- aroundSigGenes[order(aroundSigGenes)]
  #for genes having cassette exons and the -log10(P) is higher than 75
  #superSigGenes <- unique(c(rownames(which(fsmatrix >= 75, arr.ind=TRUE)), rownames(which(lsmatrix >= 75, arr.ind=TRUE))))
  #plotGenenames <- superSigGenes[order(superSigGenes)]
  
  #filename(AT1G01010)
  for(filename in plotGenenames[1:15]){
    plotcExonExp(chr, filename, fs_threshold = 5.68, ls_threshold = 5.70)
  }
}



#************************************************************* plot exp part *************************************************************#
#plot 4 env separately in 4 panel
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

plotExpEnvSep <- function(rawexp, newexp, probes_dir, exonID, uniqueExon, ind_tu, nprobes){
  #par(mfrow = c(4, 1), pty = "m", oma = c(5, 3, 5, 0.5))
  for(env in 1:4){
    ind_env <- which(as.numeric(menvironment) == env)
    
    if(env == 1) par(fig = c(0, 1, 1-0.25*env, 1.25-0.25*env), oma = c(5, 3, 5, 0.5),  mar = c(0, 4, 0, 2))
    else par(fig = c(0, 1, 1-0.25*env, 1.25-0.25*env), oma = c(5, 3, 5, 0.5),  mar = c(0, 4, 0, 2), new = TRUE)
    plot(c(0.5, nprobes + 0.5), c(min(newexp[ ,ind_env]) - 0.2, max(newexp[ ,ind_env]) + 0.2), xaxt = 'n', xlab = "", ylab = levels(menvironment)[env], cex.axis = 1, cex.lab = 1.5, col.lab = env, las = 1, mgp = c(2.25, 0.5, 0), tck = -0.017, t = "n")
    if(env == 1){
      title(main = filename, cex.main = 2.5, xlab = "Probes", mgp = c(3, 0.5, 0), cex.lab = 1.5, outer = TRUE)
      title(ylab = "Expression Intensity", mgp = c(1, 0.5, 0), cex.lab = 1.5, outer = TRUE)
    }
    
    for(p in 1:nprobes){
      #background for introns
      if(!p %in% ind_tu){
        rect((p - 0.5), -3, (p + 0.5), max(newexp[ ,ind_env])* 2, col = grey(0.85), border = "transparent")
      }
      if(p %in% probes_dir){
        points(rep(p, length(ind_env)) + runif(length(ind_env), min = -0.05, max = 0.05), newexp[p,ind_env], t = 'p', col = env, pch = 20, cex = 0.75)
      }
    }
    for(exon in uniqueExon){
      #ind <- judge which probe in exonID is of current exon name (T/F)
      ind <- rawexp[exonID, "tu"] == exon
      #cat(as.character(exon), "has probes", exonID[ind], "\n")
      if(length(exonID[ind]) >= 3) text(x = median(which(rawexp[ ,"tu"] == exon)), y = max(newexp[ ,ind_env])+0.15, labels = exon, col = "magenta4", cex = 1)
      
      #use mean/median to test cassette
      lines(c(min(which(rawexp[ ,"tu"] == exon))-0.5, max(which(rawexp[ ,"tu"] == exon))+0.5), c(mean(unlist(newexp[exonID[ind],ind_env])), mean(unlist(newexp[exonID[ind],ind_env]))), col = env, lwd = 2)       
      
      cExoninEnv <- unlist(lapply(strsplit(rownames(cematrix)[which(cematrix[ ,env] >= ce_threshold)][grepl(filename, rownames(cematrix)[which(cematrix[ ,env] >= ce_threshold)])], "_"), "[[", 2))
      if(exon %in% cExoninEnv){
        text(x = median(which(rawexp[ ,"tu"] == exon)), y = max(newexp[ ,ind_env])-0.25, labels=round(cematrix[paste0(filename, "_", exon),env], digits = 1), col = env, cex = 0.8)
      }
    }
    box()
  }
  axis(1, at = probes_dir, labels = row.names(rawexp)[probes_dir], mgp=c(2.25, 0.5, 0), cex.axis = 1, las = 2, tck = 0.02)
}


plotcExonExp <- function(chr, filename, ce_threshold){
  rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  newexp <- rawexp[ ,17:164]
  
  probes_dir <- probesDir(rawexp)
  #cat(filename, "\nprobeDir:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #cat("exons of right direction:", exonID, "\n")
  
  #uniqueExon <- all tu names of exon probes
  uniqueExon <- unique(rawexp[exonID,"tu"])
  #cat("tu names:", as.character(uniqueExon), "\n")
  
  ind_tu <- grep("tu", rawexp[ ,"tu"])
  nprobes <- nrow(rawexp)
  
  png(filename = paste0("Data/cassetteExon/chr", chr, "/", filename, "_ceExp_allind_", ce_threshold, "_4s.png"), width = 960, height = 1728, bg = "white")
  plotExpEnvSep(rawexp, newexp, probes_dir, exonID, uniqueExon, ind_tu, nprobes)
  dev.off()
}


ce_threshold = 5.86
for(chr in 1:1){
  cematrix <- read.table(paste0("Data/cassetteExon/cassetteExon_chr", chr, "_allind.txt"), row.names=1, header=T)
  
  #for all genes which have cassette exons in at least one env
  plotGenenames <- unique(unlist(lapply(strsplit(rownames(which(cematrix >= ce_threshold, arr.ind=T)), "_"), "[[", 1)))
  #for genes having cassette exons and the -log10(P) is round the ce_threshold
  #plotGenenames <- unique(unlist(lapply(strsplit(rownames(which(round(cematrix, digits=2) >= ce_threshold & round(cematrix, digits=2) < ce_threshold+0.1, arr.ind=T)), "_"), "[[", 1)))
  #for genes having cassette exons and the -log10(P) is higher than 75
  #plotGenenames <- unique(unlist(lapply(strsplit(rownames(which(cematrix > 75, arr.ind=T)), "_"), "[[", 1)))
  
  #filename(AT1G01010)
  for(filename in plotGenenames[1:15]){
    plotcExonExp(chr, filename, ce_threshold = 5.86)
  }
}
