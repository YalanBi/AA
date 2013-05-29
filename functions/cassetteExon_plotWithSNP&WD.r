#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 21-05-2013
# first written: 19-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#*************************************************************** basic part ***************************************************************#
#plot 4 env separately in 4 panel
setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]
#load names of SNP probes
SNPs <- read.table("snp ids.txt")[ ,1]
#load package scales to use function "alpha"
library("scales")


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

plotExpEnvSep <- function(filename, rawexp, newexp, probes_dir, exonID, uniqueExon, ind_tu, nprobes, ce_threshold){
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
        rect((p - 0.5), -3, (p + 0.5), max(newexp)* 2, col = grey(0.85), border = "transparent")
      }
      if(rownames(newexp)[p] %in% SNPs){
        points(rep(p, length(ind_env)) + runif(length(ind_env), min = -0.05, max = 0.05), newexp[p,ind_env], t = 'p', col = alpha(env, 0.25), pch = 20, cex = 0.5)
        text(x = p, y = min(newexp[ ,ind_env]) - 0.05, labels = "SNP", col = alpha("lemonchiffon4", 0.8), cex = 0.7)
      } else if(p %in% probes_dir){
        points(rep(p, length(ind_env)) + runif(length(ind_env), min = -0.05, max = 0.05), newexp[p,ind_env], t = 'p', col = env, pch = 20, cex = 0.75)
      } else{
        points(rep(p, length(ind_env)) + runif(length(ind_env), min = -0.05, max = 0.05), newexp[p,ind_env], t = 'p', col = alpha(env, 0.25), pch = 20, cex = 0.5)
      }
    }
    
    for(exon in uniqueExon){
      #ind <- row number of probes of current exon name, no SNP, right dir, number of rows in rawexp
      ind <- probes_dir[!rownames(newexp)[probes_dir] %in% SNPs][rawexp[probes_dir[!rownames(newexp)[probes_dir] %in% SNPs], "tu"] == exon]
      #cat(as.character(exon), "has probes", exonID[ind], "\n")
      if(length(ind) >= 3) text(x = median(which(rawexp[ ,"tu"] == exon)), y = max(newexp[ ,ind_env])+0.15, labels = exon, col = "magenta4", cex = 1)
      
      #use mean/median to test cassette
      lines(c(min(which(rawexp[ ,"tu"] == exon))-0.5, max(which(rawexp[ ,"tu"] == exon))+0.5), c(mean(unlist(newexp[ind, ind_env])), mean(unlist(newexp[ind, ind_env]))), col = env, lwd = 2)       
      
      cExoninEnv <- unlist(lapply(strsplit(rownames(cematrix)[which(cematrix[ ,env] >= ce_threshold)][grepl(filename, rownames(cematrix)[which(cematrix[ ,env] >= ce_threshold)])], "_"), "[[", 2))
      if(exon %in% cExoninEnv){
        text(x = median(which(rawexp[ ,"tu"] == exon)), y = max(newexp[ ,ind_env])-0.25, labels=round(cematrix[paste0(filename, "_", exon),env], digits = 1), col = env, cex = 0.8)
      }
    }
    box()
  }
  axis(1, at = (1:nprobes)[!1:nprobes %in% probes_dir], labels = rep("WD", nprobes-length(probes_dir)), cex.axis = 1, las = 2, tck = 0.01)
  axis(1, at = probes_dir, labels = row.names(rawexp)[probes_dir], mgp=c(2.25, 0.5, 0), cex.axis = 1, las = 2, tck = 0.02)
}


plotcExonExp <- function(chr, filename, ce_threshold){
  if(any(grepl("SNP", dir(paste0("Data/chr", chr, "_norm_hf_cor/"))[grepl(filename, dir(paste0("Data/chr", chr, "_norm_hf_cor/")))]))){
    cat("I have SNP in, let's read 'SNPin' file!\n")
    rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, "_SNPin.txt"), row.names=1, header=T)
  } else rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  newexp <- rawexp[ ,17:164]
  
  probes_dir <- probesDir(rawexp)
  #cat(filename, "\nprobeDir:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir, "tu"])]
  #cat("exons of right direction:", exonID, "\n")
  
  #uniqueExon <- all tu names of exon probes
  uniqueExon <- unique(rawexp[exonID, "tu"])
  #cat("tu names:", as.character(uniqueExon), "\n")
  
  ind_tu <- grep("tu", rawexp[ ,"tu"])
  nprobes <- nrow(rawexp)
  
  st <- proc.time()[3]
  png(filename = paste0("Data/cassetteExon/chr", chr, "/", filename, "_ceExp_allind_", ce_threshold, "_4s_SNP.png"), width = 960, height = 1728, bg = "white")
  plotExpEnvSep(filename, rawexp, newexp, probes_dir, exonID, uniqueExon, ind_tu, nprobes, ce_threshold)
  dev.off()
  et <- proc.time()[3]
  cat(filename, "finished in", et-st, "s\n")
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
