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

#cassette exon test
testCassette <- function(expdata = rawexp[exonID, ind_env + 16], ind, use){
  exonprobe <- apply(expdata[ind, ], 2, use)
  #cat(exonprobe, "\n")
  otherprobe <- apply(expdata[!ind, ], 2, use)
  #cat(otherprobe, "\n")
  return(-log10(t.test(exonprobe, otherprobe, alternative="less")$p.value))
}


#*************************************************************** load part ***************************************************************#
#load exp genes
load(file="Data/fullModeMapping/expGenes.Rdata")


#*************************************************************** test part ***************************************************************#
gene_count <- c(0, 0, 0, 0, 0)
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
    
    #at least 3 exons in a gene!!!
    if(length(uniqueExon) >= 3){
      #cat("we have >= 3 exons!\n")
      gene_count[chr] <- gene_count[chr] + 1
      for(exon in uniqueExon[2: (length(uniqueExon)-1)]){
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
            res <- c(res, testCassette(expdata = rawexp[exonID, ind_env + 16], ind, use = median)) #********** change!!! **********#
            
            #t.test once, counter plus 1!!! to calculate the 
            cnt <- cnt + 1
          }
       resmatrix <- rbind(resmatrix, res)
        }
      }
    }
  }
  cnt_mem <- c(cnt_mem, cnt)
  
  rownames(resmatrix) <- rownameList
  colnames(resmatrix) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
  write.table(resmatrix, file=paste0("Data/cassetteExon/cassetteExon_chr", chr, "_median_D2p.txt"), sep="\t") #change!!! 
  
  et <- proc.time()[3]
  cat("and finished in", et-st, "s\n")
}

sum(cnt_mem)
#93820
gene_count#(for median)
[1] 2395 1365 1836 1406 2165
#sum(gene_count)=9167

#count ngene having

#************************************************************* summary part *************************************************************#
#mean or median, change!!!
res <- NULL
for(chr in 1:5){
  cematrix <- read.table(paste0("Data/cassetteExon/cassetteExon_chr", chr, "_mean_D2p.txt"), row.names=1, header=TRUE) #********** change!!! **********#
  rr <- NULL
  for(env in 1:4){
    rr <- c(rr, length(which(cematrix[ ,env] >= -log10(0.05/sum(cnt_mem))))/sum(cnt_mem/4))
  }
  res <- rbind(res, rr)
}
res <- rbind(res, colSums(res))
rownames(res) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "Sum")
colnames(res) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
write.table(res, file=paste0("Data/cassetteExon/cassetteExon_sum_mean_D2p.txt"), sep="\t") #********** change!!! **********#










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
      rect((p-0.5), -3, (p+0.5), max(newexp)* 2.1, col=grey(0.85), border = "transparent")
    }
    if(p %in% probes_dir){
      points(rep(p,148)+0.06*as.numeric(menvironment)-0.15, newexp[p,], t='p', col=menvironment, pch=20, cex=0.7)
    }
  }
  for(exon in uniqueExon){
    #ind <- judge which probe in exonID is of current exon name (T/F)
    ind <- rawexp[exonID, "tu"] == exon
    #cat(as.character(exon), "has probes", exonID[ind], "\n")
    
    for(env in 1:4){
      ind_env <- which(as.numeric(menvironment) == env)
      
      #use mean/median to test cassette
      lines(c(min(which(rawexp[,"tu"] == exon))-0.5+0.1*env, max(which(rawexp[,"tu"] == exon))+0.1*env), c(mean(unlist(newexp[exonID[ind], ind_env])), mean(unlist(newexp[exonID[ind], ind_env]))), col=env, lwd=2)       
      
      cExoninEnv <- unlist(lapply(strsplit(row.names(cematrix[which(cematrix[,env]>=ce_threshold),])[grepl(filename, row.names(cematrix[which(cematrix[,env]>=ce_threshold),]))],"_"),"[[",2))
      if(exon %in% cExoninEnv){
        text(x=median(exonID[ind]), y=max(newexp)-0.1*env, labels=exon, col=env, cex=1)
      }
    }
  }
  axis(1, at=probes_dir, labels=row.names(rawexp)[probes_dir], cex.axis=1, las=2, tck=0.005)
  box()
  dev.off()
}


use = "median"
for(chr in 1:1){
  cematrix <- read.table(paste0("Data/cassetteExon/cassetteExon_chr1_", use, "_D2p.txt"), row.names=1, header=T)
  
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
