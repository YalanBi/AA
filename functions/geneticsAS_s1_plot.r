#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 05-06-2013
# first written: 03-06-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#************************************************************* plot exp part *************************************************************#
#plot 4 env separately in 4 panel
setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]
#load genotype file
geno <- read.table("refined map/genotypes.txt",sep="\t", row.names=1, header=TRUE)



ce_threshold = 5.86; threshold_qtl = 8.0; threshold_int = 11.6; cutoffnProbe = 2; #cutoffratio = 0.6#(cutoffratio>=0.6; cutoffnProbe >= 2 probes)
#*************************************************************** load part ***************************************************************#
#load main/int QTL genes
load(file=paste0("Data/geneticsAS/genelist_QTL", threshold_qtl, "_Int",threshold_int, "_np", cutoffnProbe, ".Rdata"))


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


plotExpEnvSep <- function(chr, filename, rawexp, newexp, probes_dir, exonID, uniqueExon, ind_tu, nprobes, ASexon, genoProbes, markerToDraw, int, ce_threshold){
  #par(mfrow = c(4, 1), pty = "m", oma = c(5, 3, 5, 0.5))
  for(env in 1:4){
    ind_env <- which(as.numeric(menvironment) == env)
    
    if(env == 1) par(fig = c(0, 1, 1-0.25*env, 1.25-0.25*env), oma = c(5, 3, 5, 0.5),  mar = c(0, 4, 0, 2))
    else par(fig = c(0, 1, 1-0.25*env, 1.25-0.25*env), oma = c(5, 3, 5, 0.5),  mar = c(0, 4, 0, 2), new = TRUE)
    plot(c(0.5, nprobes + 0.5), c(min(newexp) - 0.2, max(newexp) + 0.2), xaxt = 'n', xlab = "", ylab = levels(menvironment)[env], cex.axis = 1, cex.lab = 1.5, col.lab = env, las = 1, mgp = c(2.25, 0.5, 0), tck = -0.017, t = "n")
    if(env == 1){
      title(main = filename, cex.main = 2.5, xlab = "Probes", mgp = c(3, 0.5, 0), cex.lab = 1.5, outer = TRUE)
      title(ylab = "Expression Intensity", mgp = c(1, 0.5, 0), cex.lab = 1.5, outer = TRUE)
    }
    
    for(p in 1:nprobes){
      #background for introns
      if(!p %in% ind_tu) rect((p - 0.5), -3, (p + 0.5), max(newexp[ ,ind_env])* 2, col = grey(0.85), border = "transparent")
      
      if(p %in% genoProbes){
        if(length(grep(filename, genesIAS[[paste0("chr", chr)]][[env]], value=T)) > 0) points(rep(p, length(ind_env)) + runif(length(ind_env), min = -0.05, max = 0.05), newexp[p,ind_env], t = 'p', col = geno[ ,markerToDraw[1]]+4, pch = 20, cex = 0.75)
        else points(rep(p, length(ind_env)) + runif(length(ind_env), min = -0.05, max = 0.05), newexp[p,ind_env], t = 'p', col = env, pch = 20, cex = 0.75)
        points(p, mean(unlist(newexp[p,ind_env[geno[ind_env,markerToDraw[1]] == 1]])), pch="*", cex=2, col=5)
        points(p, mean(unlist(newexp[p,ind_env[geno[ind_env,markerToDraw[1]] == 2]])), pch="*", cex=2, col=6)
        
        if(int[p,markerToDraw[1]] >= threshold_int && env %in% grep(paste0(filename, "_", rawexp[p, "tu"]), genesIAS[[paste0("chr", chr)]], value=F)){
          text(x = p, y = min(newexp)+0.15, labels = round(int[p,markerToDraw[1]], digits = 1), col = env, cex = 0.8)
        }
      } else if(p %in% probes_dir) points(rep(p, length(ind_env)) + runif(length(ind_env), min = -0.05, max = 0.05), newexp[p,ind_env], t = 'p', col = env, pch = 20, cex = 0.75)
    }
        
    for(exon in uniqueExon){
      #ind <- judge which probe in exonID is of current exon name (T/F)
      ind <- rawexp[exonID, "tu"] == exon
      #cat(as.character(exon), "has probes", exonID[ind], "\n")
      #if(length(exonID[ind]) >= 3) text(x = median(which(rawexp[ ,"tu"] == exon)), y = max(newexp)-0.1, labels = exon, col = "magenta4", cex = 1)
      
      #use mean/median to test cassette
      lines(c(min(which(rawexp[ ,"tu"] == exon))-0.5, max(which(rawexp[ ,"tu"] == exon))+0.5), c(mean(unlist(newexp[exonID[ind],ind_env])), mean(unlist(newexp[exonID[ind],ind_env]))), col = env, lwd = 2)       
      
      cExoninEnv <- unlist(lapply(strsplit(rownames(cematrix)[which(cematrix[ ,env] >= ce_threshold)][grepl(filename, rownames(cematrix)[which(cematrix[ ,env] >= ce_threshold)])], "_"), "[[", 2))
      if(exon %in% cExoninEnv){
        text(x = median(which(rawexp[ ,"tu"] == exon)), y = max(newexp)-0.25, labels=paste0("ce=", round(cematrix[paste0(filename, "_", exon),env], digits = 1)), col = env, cex = 0.8)
      }
      if(exon %in% ASexon){
        text(x = median(which(rawexp[ ,"tu"] == exon)), y = max(newexp)-0.4, labels=paste0("m=", markerToDraw[1]), col = "magenta4", cex = 0.8)
      }
    }
    box()
  }
  axis(1, at = probes_dir, labels = row.names(rawexp)[probes_dir], mgp=c(2.25, 0.5, 0), cex.axis = 1, las = 2, tck = 0.02)
}


plotcExonExp <- function(chr, filename, ce_threshold){
  rawexp <- read.table(paste0("Data/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  newexp <- rawexp[ ,17:164]
  int <- read.table(paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/", filename, "_FM_Int.txt"), row.names=1, header=T)
  #cat(" int loading succeed!\n")
  
  probes_dir <- probesDir(rawexp)
  #cat(filename, "\nprobeDir:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #cat("exons of right direction:", exonID, "\n")
  
  #uniqueExon <- all tu names of exon probes
  uniqueExon <- unique(rawexp[exonID,"tu"])
  #cat("tu names:", as.character(uniqueExon), "\n")
  
  ind_tu <- grep("tu", rawexp[ ,"tu"])
  nprobes <- nrow(rawexp)
  
  genoProbes <- NULL
  for(env in 1:4){
    if(length(grep(filename, genesIAS[[paste0("chr", chr)]][[env]], value=T)) > 0 && length(genoProbes) == 0){
      ASexon <- unique(unlist(lapply(strsplit(grep(filename, genesIAS[[paste0("chr", chr)]][[env]], value=T), "_"), "[[", 2)))
      genoProbes <- exonID[rawexp[exonID,"tu"] == ASexon]
      consSigIntMarkers <- as.numeric(unlist(lapply(strsplit(grep(filename, genesIAS[[paste0("chr", chr)]][[env]], value=T), "_"), "[[", 3)))
      markerToDraw <- consSigIntMarkers[which(int[genoProbes,consSigIntMarkers] == max(int[genoProbes,consSigIntMarkers]), arr.ind=T)[,2]]
      #cat("env", env, ": we are markers, at which there're >= 2 probes with sig Int,", markerToDraw, "\n")
    }
  }
  
  png(filename = paste0("Data/geneticsAS/plotGeneticsAS/", filename, "_QTL", threshold_qtl, "_Int", threshold_int, "_np", cutoffnProbe, "_4s.png"), width = 960, height = 1728, bg = "white")
  plotExpEnvSep(chr, filename, rawexp, newexp, probes_dir, exonID, uniqueExon, ind_tu, nprobes, ASexon, genoProbes, markerToDraw, int, ce_threshold)
  dev.off()
}


for(chr in 1:5){
  cematrix <- read.table(paste0("Data/cassetteExon/cassetteExon_chr", chr, "_allind.txt"), row.names=1, header=T)
  
  #for all genes which have cassette exons in at least one env
  if(length(unlist(genesIAS[[paste0("chr", chr)]])) > 0){
    plotGenenames <- sort(unique(unlist(lapply(strsplit(unlist(genesIAS[[paste0("chr", chr)]]), "_"), "[[", 1))))
    
    #filename(AT1G01010)
    for(filename in plotGenenames){
      plotcExonExp(chr, filename, ce_threshold = 5.86)
    }
  }
}
