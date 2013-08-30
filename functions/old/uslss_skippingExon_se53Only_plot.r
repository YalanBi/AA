
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
