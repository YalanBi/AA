#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 19-06-2013
# first written: 21-05-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#********************************************** This is used to plot alternative splicing at 5' **********************************************#
#plot 4 env separately in 4 panel
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

plotExpEnvSep <- function(filename, rawexp, newexp, probes_dir, exonID, uniqueExon, ind_tu, nprobes, testInd, thr = ps_threshold){
  #par(mfrow = c(4, 1), pty = "m", oma = c(5, 3, 5, 0.5))
  
  #which tu is in rownames(test result matrix)
  testTuID <- unlist(lapply(strsplit(rownames(psmatrix)[grepl(filename, rownames(psmatrix))], "_"), "[[", 2))
  
  for(env in 1:4){
    ind_env <- which(as.numeric(menvironment) == env)
    
    if(env == 1) par(fig = c(0, 1, 1-0.2*env, 1.2-0.2*env), oma = c(5, 3, 5, 0.5),  mar = c(0, 4, 0, 2))
    else par(fig = c(0, 1, 1-0.25*env, 1.25-0.25*env), oma = c(5, 3, 5, 0.5),  mar = c(0, 4, 0, 2), new = TRUE)
    plot(c(0.5, nprobes + 0.5), c(min(newexp[ ,ind_env]) - 0.2, max(newexp[ ,ind_env]) + 0.2), xaxt = 'n', xlab = "", ylab = levels(menvironment)[env], cex.axis = 1, cex.lab = 1.5, col.lab = env, las = 1, mgp = c(2.25, 0.5, 0), tck = -0.017, t = "n")
    if(env == 1){
      title(main = paste0(filename, "_5'AS"), cex.main = 2.5, xlab = "Probes", mgp = c(3, 0.5, 0), cex.lab = 1.5, outer = TRUE)
      title(ylab = "Expression Intensity", mgp = c(1, 0.5, 0), cex.lab = 1.5, outer = TRUE)
      
      axis(3, at = mean(which(rawexp[ ,"tu"] == tu)), labels = tu, mgp=c(2.25, 0.5, 0), tick = FALSE, line = 0)
    }
    
    if(any(testTuID == 1)){
      
      rect(0.5, -3, (exonID[ind][which.max(dff)] + exonID[ind][which.max(dff) + 1]) * 0.5, 20, density = 15, angle = 30, col = "azure2", border = "azure3")
      rect((exonID[ind][which.max(dff)] + exonID[ind][which.max(dff) + 1]) * 0.5, -3, max(which(rawexp[, "tu"] == uniqueExon[1])) + 0.5, 20, density = 15, angle = 150, col = "azure2", border = "azure3")
      
      #use mean/median to test AS in 5'|3' site
      lines(c(0.5, (exonID[ind][which.max(dff)] + exonID[ind][which.max(dff) + 1]) * 0.5), c(mean(unlist(newexp[exonID[ind][1:which.max(dff)],ind_env])), mean(unlist(newexp[exonID[ind][1:which.max(dff)],ind_env]))), col = env, lwd = 2)
      lines(c((exonID[ind][which.max(dff)] + exonID[ind][which.max(dff) + 1]) * 0.5, max(which(rawexp[, "tu"] == uniqueExon[1]))+0.5), c(mean(unlist(newexp[exonID[ind][-(1:which.max(dff))],ind_env])), mean(unlist(newexp[exonID[ind][-(1:which.max(dff))],ind_env]))), col = env, lwd = 2)
      
      if(psmatrix[rownames(psmatrix)[grepl(filename, rownames(psmatrix))][nc],env] >= thr){
        text(x = median(which(rawexp[ ,"tu"] == uniqueExon[1])), y = max(newexp[ ,ind_env])+0.15, labels = round(psmatrix[rownames(psmatrix)[grepl(filename, rownames(psmatrix))][nc],env], digits = 1), col = env, cex = 0.8)
      }
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
    
    for(e in 1:length(uniqueExon)){
      if(! e %in% testTuID){
        ind <- rawexp[exonID, "tu"] == uniqueExon[e]
        lines(c(min(which(rawexp[,"tu"] == uniqueExon[e])) - 0.5, max(which(rawexp[,"tu"] == uniqueExon[e])) + 0.5), c(median(unlist(newexp[exonID[ind],ind_env])), median(unlist(newexp[exonID[ind],ind_env]))), col = env, lwd = 2)
      }
    }
    
    box()
  }
  axis(1, at = probes_dir, labels = row.names(rawexp)[probes_dir], mgp=c(2.25, 0.5, 0), cex.axis = 1, las = 2, tck = 0.02)
}


plotcExonExp <- function(chr, filename, ps_threshold){
  rawexp <- read.table(paste0("Data/Raw/chr", chr, "_norm_hf_cor/", filename, ".txt"), row.names=1, header=T)
  newexp <- rawexp[ ,17:164]
  
  probes_dir <- probesDir(rawexp)
  #cat(filename, "\nprobeDir:", probes_dir, "\n")
  exonID <- probes_dir[grepl("tu", rawexp[probes_dir,"tu"])]
  #cat("exons of right direction:", exonID, "\n")
  
  #uniqueExon <- all tu names of exon probes
  uniqueExon <- unique(grep("tu", rawexp[ ,"tu"], value=TRUE))
  #cat("tu names:", as.character(uniqueExon), "\n")
  
  ind_tu <- grep("tu", rawexp[ ,"tu"])
  nprobes <- nrow(rawexp)
  
  testInd <- grep(uniqueExon[1], rawexp[exonID,"tu"])
  dff <- NULL
  for(n in 2:length(ind)){
    dff <- c(dff, sum(abs(rawexp[exonID[ind][n], 17:164] - rawexp[exonID[ind][n-1], 17:164])))
    #cat("  difference between p", exonID[ind][n], "and p", exonID[ind][n-1], "is", sum(rawexp[exonID[ind][n], 17:164]) - sum(rawexp[exonID[ind][n-1], 17:164]), "\n")
  }
  
  png(filename = paste0("Data/53terminalAS/plot_53ps/", filename, "_wtest_", ps_threshold, "_4s.png"), width = 960, height = 1728, bg = "white")
  plotExpEnvSep(filename, rawexp, newexp, probes_dir, exonID, uniqueExon, ind_tu, nprobes, thr = ps_threshold)
  dev.off()
}


#ps -> partly splicing
as3Thre <- 5.09
for(chr in 1:5){
  load(file = paste0("Data/53terminalAS/3'AS_chr", chr, "_genenameList.Rdata"))
  genenames <- as5ResList[[5]]
  
  psmatrix <- read.table(paste0("Data/53terminalAS/53terminalAS_chr", chr, "_ttest.txt"), row.names=1, header=T)
  
  #for all genes which have cassette exons in at least one env
  plotGenenames <- sort(unique(unlist(lapply(strsplit(rownames(which(psmatrix >= ps_threshold, arr.ind=T)), "_"), "[[", 1))))
  #for genes having cassette exons and the -log10(P) is round the ps_threshold
  #plotGenenames <- sort(unique(unlist(lapply(strsplit(rownames(which(round(psmatrix, digits=2) >= ps_threshold & round(psmatrix, digits=2) < ps_threshold+0.1, arr.ind=T)), "_"), "[[", 1))))
  #for genes having cassette exons and the -log10(P) is higher than 75
  #plotGenenames <- sort(unique(unlist(lapply(strsplit(rownames(which(psmatrix > 75, arr.ind=T)), "_"), "[[", 1))))
  
  #filename(AT1G01010)
  for(filename in plotGenenames){
    plotcExonExp(chr, filename, ps_threshold = 5.03)
  }
}



#*********************************************************** find best examples ***********************************************************#
getOne <- function(aa, cond = 1, cutoff=30){
  bb <- aa[which(aa[ ,cond] > cutoff), ]
  l10 <- sort(apply(bb[ ,(1:4)[-cond]], 1, sum), index.return=T)$ix[1:10]#sort(apply(bb[ ,-cond], 1, sum), index.return=T)$ix[1:10]
  bb[l10, ]
}

getOne <- function(aa, cond = 1, cutoff=30){
  bb <- aa[which(aa[ ,cond] > cutoff), ]
  l10 <- NULL
  for(r in 1:nrow(bb)){
    if(any(bb[r, -cond] <= (cutoff - 15))) l10 <- c(l10, r)
  }
  bb[l10, ]
}

for(chr in 1:5){
  #psmatrix <- read.table(paste0("Data/53terminalAS/53terminalAS_chr", chr, "_ttest.txt"), row.names=1, header=T)
  psmatrix <- read.table(paste0("Data/53terminalAS/53terminalAS_chr", chr, "_wtest.txt"), row.names=1, header=T)
  
  posPerfectEg <- NULL
  for(env in 1:4){
    posPerfectEg <- c(posPerfectEg, unlist(lapply(strsplit(rownames(getOne(psmatrix, env, 30)), "_"), "[[", 1)))
  }
  cat("I'm chr", chr, ", I have\n", posPerfectEg, "\n")
  for(filename in sort(unique(posPerfectEg))){
    #chr <- gsub("AT", "", unlist(lapply(strsplit(filename, "G"), "[[", 1)))
    plotcExonExp(chr, filename, ps_threshold = 5.03)
  }
}
