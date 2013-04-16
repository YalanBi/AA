#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 11-04-2013
# first written: 22-03-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")
#load environment file
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]
#load genomap
ann_m <- read.table("refined map/map.txt")

############################################################## probes selection ##############################################################
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





#*************************************************************** counting part ***************************************************************#
#test for intron retention: TRUE / FALSE
checkRetentionTF <- function(rawexp, newexp, intronID, exonExpID, nIn, ind_env, threshold){
  
  #retentionList <- a list of retained intron names of current gene under this environment
  retained <- FALSE
  for(i in unique(rawexp[intronID,"tu"])){
    
    #probes4i <- probes that of current intron names and of the right direction
    probes4i <- intronID[which(rawexp[intronID,"tu"]==i)]
    #cat(i, ":", probes4i, "\n")
    
    #to decide the min number of intron probes
    if(length(probes4i) >= nIn){
      
      #compare mean(current intron probes) and mean(expressed exons which are in the right direction)
      if(mean(unlist(newexp[probes4i, ind_env])) >= mean(unlist(newexp[exonExpID, ind_env]))){
        retained <- TRUE
        #cat("means of", i, "is higher than expressed exons\n")
      }else{ 
        #t.test (current intron probes) against (expressed exons which are in the right direction), we want unsignificant p-value
        if(-log10(t.test(newexp[probes4i, ind_env], newexp[exonExpID, ind_env])$p.value) <= threshold){
          retained <- TRUE
          #cat("means of", i, "is lower than expressed exons, but not sig different.\n")
        }
      }
    }
  }
  return(retained)
}

################################################################# count by chr #################################################################
countReBYchr <- function(chr, menvironment, nIn, cutoff, threshold){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")

  #number of genes that expressed in each environmet and have at least 1 intron probe in the right direction!!!
  countGene <- c(0,0,0,0)
  #number of genes that has retention in each environmet
  countRetention <- c(0,0,0,0)
  
  #only find genes containing more than 4 probes in total
  for(filename in gsub("_QTL", "", dir(location)[which(grepl("_QTL", dir(location)))])){
    rawexp <- read.table(paste0(location, filename), row.names=1, header=T)
    newexp <- rawexp[,17:164]
    
    probes_dir <- probesDir(rawexp)
    cat(filename, "\nprobeDir:", probes_dir, "\n")
    #exonID <- exons of right direction
    exonID <- probes_dir[which(grepl("tu", rawexp[probes_dir, "tu"]))]
    cat("exons:", exonID, "\n")
    #intronID <- introns of right direction
    intronID <- probes_dir[which(grepl("intron", rawexp[probes_dir, "tu"]))]
    cat("introns:", intronID, "\n")
    
    #there must be at least 1 intron probe, then we test for intron retention
    if(length(intronID) > 0){
      for(env in 1:4){
        ind_env <- which(as.numeric(menvironment) == env)
        retainedORnot <- FALSE
        
        #check whether current gene expressed or not under this environment, True: continue test intron retention
        #cutoff <- cutoff used to check whether mean(all exons of right direction) is high enough to be regatded as expressed or not
        if(length(exonID) > 0 && mean(unlist(newexp[exonID, ind_env])) >= cutoff){
          cat("means of exons in env", env, "higher than cutoff=", cutoff, "\n")
          
          exonExpID <- chooseExpExon(exonID, newexp, ind_env, cutoff)
          cat("expressed exons:", exonExpID, "\n")
          
          retainedORnot <- checkRetentionTF(rawexp, newexp, intronID, exonExpID, nIn, ind_env, threshold)
          cat("introns retained in env", env, retainedORnot, "\n")
          
          countGene[env] <- countGene[env] + 1
        }
        else{
          cat("means of exons in env", env, "lower than cutoff=", cutoff, "\n")
        }
        
        if(retainedORnot){
          countRetention[env] <- countRetention[env] + 1
        }
      }
    } else{
      cat("no intron!\n")
    }
    
  }
  res <- rbind(countRetention, countGene, countRetention/countGene)
  rownames(res) <- c(paste0("chr", chr, "_nRetained"), paste0("chr", chr, "_nExpGene"), paste0("chr", chr, "_ratio"))
  return(res)
}


resmatrix <- NULL
for(chr in 1:5){
  resmatrix <- rbind(resmatrix, countReBYchr(chr, menvironment, nIn = 1, cutoff = 4.5, threshold = 10))
  colnames(resmatrix) <- c("6H", "Dry_AR", "Dry_Fresh", "RP")
}
resmatrix
cat("nIntronProbes1; expCutoff4.5; pThres10\n")

resmatrix <- rbind(resmatrix,resmatrix[1,]+resmatrix[4,]+resmatrix[7,]+resmatrix[10,]+resmatrix[13,], resmatrix[2,]+resmatrix[5,]+resmatrix[8,]+resmatrix[11,]+resmatrix[14,], (resmatrix[1,]+resmatrix[4,]+resmatrix[7,]+resmatrix[10,]+resmatrix[13,])/(resmatrix[2,]+resmatrix[5,]+resmatrix[8,]+resmatrix[11,]+resmatrix[14,]))
rownames(resmatrix)[16:18] <- c("genome_nRetained", "genome_nExpGene", "genome_ratio")
write.table(resmatrix, file="Data/intronRetention/intronRetention_ratioSum_nIn1_Exp4.5_reThres10.txt")



#*************************************************************** plot part ***************************************************************#
#test for intron retention: a list of introns
checkRetentionList <- function(rawexp, newexp, intronID, exonExpID, nIn, ind_env, threshold){
  
  #retentionList <- a list of retained intron names of current gene under this environment
  retentionList <- NULL
  for(i in unique(rawexp[intronID,"tu"])){
    
    #probes4i <- probes that of current intron names and of the right direction
    probes4i <- intronID[which(rawexp[intronID,"tu"]==i)]
    #cat(i, ":", probes4i, "\n")
    
    #to decide the min number of intron probes
    if(length(probes4i) >= nIn){
      
      #compare mean(current intron probes) and mean(expressed exons which are in the right direction)
      if(mean(unlist(newexp[probes4i, ind_env])) >= mean(unlist(newexp[exonExpID, ind_env]))){
        retentionList <- c(retentionList, i)
        #cat("means of", i, "is higher than expressed exons\n")
      } else{
        #t.test (current intron probes) against (expressed exons which are in the right direction), we want unsignificant p-value
        if(-log10(t.test(newexp[probes4i, ind_env], newexp[exonExpID, ind_env])$p.value) <= threshold){
          retentionList <- c(retentionList, i)
          #cat("means of", i, "is lower than expressed exons, but not sig different.\n")
        }
      }
    }
  }
  return(retentionList)
}

################################################################# panel 1 #################################################################
#plot for expression intensity
makePlot_Exp <- function(filename, rawexp, newexp, nprobes, intronID, exonExpID, menvironment, ind_tu, res){
  plot(c(0.5, nprobes+0.5), c(round(min(newexp)-0.5), round(max(newexp)+0.5)), xaxt='n', xlab="", ylab="Intensity", cex.axis=1, cex.lab=1.5, cex.main=2, las=1, mgp=c(3,1,0), tck=-0.017, t="n", main=filename)
  
  for(p in 1:nprobes){
    
    #background for introns
    if(!p %in% ind_tu){
      rect((p-0.5), -3, (p+0.5), max(newexp)* 2.1, col=grey(0.85), border = "transparent")
    }
    
    #which probes are showing in the panel 1: expressed exons and introns
    probesID <- unique(c(exonExpID, intronID))[order(unique(c(exonExpID, intronID)))]
    #points of expression of probes used to check intron retention
    if(p %in% probesID){
      points(rep(p,148)+0.06*as.numeric(menvironment)-0.15, newexp[p,], t='p', col=menvironment, pch=20, cex=0.7)
    }
    
    #show intron retention
    for(env in 1:4){
      #if in current environment, there is retention, show retention
      if(length(res[[env]]) > 0){
        #present means of expressed exon in current environment
        lines(c(0+0.2*env, nprobes+0.2*env), c(mean(unlist(newexp[exonExpID, which(as.numeric(menvironment) == env)])), mean(unlist(newexp[exonExpID, which(as.numeric(menvironment) == env)]))), col=as.factor(levels(menvironment))[env], lty=8, lwd=1)
        
        for(i in res[[env]]){
          ind_env <- which(as.numeric(menvironment) == env)
          probes4i <- intronID[which(rawexp[intronID,"tu"]==i)]
          lines(c(min(which(rawexp[ ,"tu"] == i)-0.5), max(which(rawexp[ ,"tu"] == i)+0.5)), c(mean(unlist(newexp[probes4i, ind_env])), mean(unlist(newexp[probes4i, ind_env]))), col=as.factor(levels(menvironment))[env], lwd=2)
        }
      }
    }
  }
  
  box();
}

################################################################# panel 2 #################################################################
#separate markers by chr
getProbesOnChr <- function(ann_m, chrs = 1){
  which(ann_m[,1] == chrs)
}

#plot for eQTL
makePlot_eQTL <- function(qtl, nprobes, ind_tu, lodThreshold = 5, lchr = chr){
  plot(c(0.5, nprobes+0.5),c(0.5, 5.5), xaxt='n', xlab="", yaxt='n', ylab="eQTL", cex.axis=1, cex.lab=1.5, las=1, mgp=c(3,1,0), t="n")
  axis(2, at=1:5, labels=paste("chr", c("I", "II", "III", "IV", "V"),sep=""), cex.axis=1, las=2, tck=-0.035, mgp=c(2,0.5,0))
  
  #background for introns
  for(p in 1:nprobes){
    if(!p %in% ind_tu){
      rect((p-0.5), -3, (p+0.5), 8, col=grey(0.85), border = "transparent")
    }
  }
  
  #line out which chr current gene are on
  abline(h=lchr, col="burlywood3", lty=8, lwd=3)
  points(-1, lchr, pch=">", cex=5, col="burlywood3")
  
  #rectangles for sig eQTL >= lodThreshold
  for(p in 1:nprobes){
    for(chrs in 1:5){
      if(any(qtl[p, getProbesOnChr(ann_m, chrs)] >= lodThreshold)){
        #cat("->",(round(max(qtl[p, getProbesOnChr(ann_m, chrs)]),d=0)+1) - lodThreshold,"\n")
        #points(p, chrs, pch=15,cex=1)
        rect((p-0.5),(chrs-0.4),(p+0.5),(chrs+0.4), col='tan2', border = "sienna3")
      }
    }
  }
  
  box();
}

################################################################# panel 3 #################################################################
#means of individuals under each environment for every probe
envTtest2 <- function(newexp, menvironment, p){
  env1 <- mean(as.numeric(newexp[p, which(as.numeric(menvironment)==1)]))
  env2 <- mean(as.numeric(newexp[p, which(as.numeric(menvironment)==2)]))
  env3 <- mean(as.numeric(newexp[p, which(as.numeric(menvironment)==3)]))
  env4 <- mean(as.numeric(newexp[p, which(as.numeric(menvironment)==4)]))
  #lgp <- -log10(c(t.test(env1, mu=mean(c(env2,env3,env4)))$p.value,t.test(env2, mu=mean(c(env1,env3,env4)))$p.value,t.test(env3, mu=mean(c(env1,env2,env4)))$p.value,t.test(env4, mu=mean(c(env1,env2,env3)))$p.value))
  #env_mean <- mean(c(env1,env2,env3,env4))
  #lgp <- -log10(c(t.test(env1, mu=env_mean)$p.value,t.test(env2, mu=env_mean)$p.value,t.test(env3, mu=env_mean)$p.value,t.test(env4, mu=env_mean)$p.value))
  return(c(env1,env2,env3,env4))
}

#load color package and design color pattern
library("colorspace", lib.loc="C:/R/win-library/2.15")
cols <- diverge_hcl(31, h=c(195,330), c = 95, l = c(20, 90), power = 1.25)

#plot of environment comparison
makePlot_Env <- function(newexp, nprobes, menvironment){
  plot(c(0.5, nprobes+0.5),c(0.5,4.5), xaxt='n', xlab="", yaxt='n', ylab="Env", cex.axis=1, cex.lab=1.5, las=1, mgp=c(3,1,0), t="n")
  axis(1, at=1:nprobes, labels=row.names(newexp), cex.axis=1, las=2, tck=0.035, mgp=c(2,1,0))
  axis(2, at=1, labels="6H", cex.axis=1, col.axis='black', las=2, tck=-0.035, mgp=c(2,0.5,0))
  axis(2, at=2, labels="DA", cex.axis=1, col.axis='red', las=2, tck=-0.035, mgp=c(2,0.5,0))
  axis(2, at=3, labels="DF", cex.axis=1, col.axis='green', las=2, tck=-0.035, mgp=c(2,0.5,0))
  axis(2, at=4, labels="RP", cex.axis=1, col.axis='blue', las=2, tck=-0.035, mgp=c(2,0.5,0))
  for(p in 1:nprobes){
    #pForCol <- round(envTtest(newexp,env,p)-0.5)+1
    pForCol <- ((envTtest2(newexp, menvironment, p) - mean(as.numeric(newexp[p,])))*7.5) + 16
    pForCol[pForCol > 31] <- 31
    pForCol[pForCol < 1] <- 1
    #cat(pForCol,"\n")
    rect((p-0.5),0.6,(p+0.5),1.4,col=cols[pForCol[1]],border = "transparent")
    rect((p-0.5),1.6,(p+0.5),2.4,col=cols[pForCol[2]],border = "transparent")
    rect((p-0.5),2.6,(p+0.5),3.4,col=cols[pForCol[3]],border = "transparent")
    rect((p-0.5),3.6,(p+0.5),4.4,col=cols[pForCol[4]],border = "transparent")
  }
  box();
}

################################################################# make 3 in 1 #################################################################
singlePlot <- function(chr, rawexp, newexp, qtl, filename, intronID, exonExpID, menvironment, ind_tu, res, lodThreshold){
  ind_tu <- grep("tu", rawexp[,"tu"])
  nprobes <- nrow(qtl)
  st <- proc.time()
  par(fig=c(0,1,0.5,1), mar=c(0,5,5,2), oma=c(0,0,0,0.5))
  makePlot_Exp(filename, rawexp, newexp, nprobes, intronID, exonExpID, menvironment, ind_tu, res)
  par(fig=c(0,1,0.3,0.5), mar=c(0,5,0,2), oma=c(0,0,0,0.5), new=T)
  makePlot_eQTL(qtl, nprobes, ind_tu, lodThreshold = 5, lchr = chr)
  par(fig=c(0,1,0,0.3), mar=c(5,5,0,2), oma=c(0,0,0,0.5), new=T)
  makePlot_Env(newexp, nprobes, menvironment)
  et <- proc.time()
  cat("Done with", filename, "after:",(et-st)[3],"secs\n")
}

######################################################### plot 20 pictures each chr #########################################################
plotReBYchr <- function(chr, menvironment, nIn, cutoff = 5, threshold = 7, lodThreshold = 4){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  nPlot <- 0
  
  #only find genes containing more than 4 probes in total
  #filename(AT1G01010)
  for(filename in gsub("_QTL.txt", "", dir(location)[which(grepl("_QTL", dir(location)))])){
    if(nPlot < 20){
      rawexp <- read.table(paste0(location, filename, ".txt"), row.names=1, header=T)
      newexp <- rawexp[,17:164]
      
      probes_dir <- probesDir(rawexp)
      #cat(filename, "\nprobeDir:", probes_dir, "\n")
      #exonID <- exons of right direction
      exonID <- probes_dir[which(grepl("tu", rawexp[probes_dir, "tu"]))]
      #cat("exons:", exonID, "\n")
      #intronID <- introns of right direction
      intronID <- probes_dir[which(grepl("intron", rawexp[probes_dir, "tu"]))]
      #cat("introns:", intronID, "\n")
      
      #store result by environment
      res <- vector("list", 4)
      
      #decide whether to make plot or not
      draw <- FALSE
      
      if(length(intronID) > 0){
        for(env in 1:4){
          ind_env <- which(as.numeric(menvironment) == env)
          retainedInList <- NULL
          
          #check whether current gene expressed or not under this environment, True: continue test intron retention
          #cutoff <- cutoff used to check whether mean(all exons of right direction) is high enough to be regatded as expressed or not
          if(length(exonID) > 0 && mean(unlist(newexp[exonID, ind_env])) >= cutoff){
            #cat("means of exons in env", env, "higher than cutoff=", cutoff, "\n")
            
            exonExpID <- chooseExpExon(exonID, newexp, ind_env, cutoff)
            #cat("expressed exons:", exonExpID, "\n")
            
            retainedInList <- checkRetentionList(rawexp, newexp, intronID, exonExpID, nIn, ind_env, threshold)
            #cat("introns in env", env, "lower than threshold=", threshold, ":", retainedInList, "\n")
            
          }
          else{
            #cat("means of exons in env", env, "lower than cutoff=", cutoff, "\n")
          }
          
          if(length(retainedInList) != 0){
            draw <- TRUE
            res[[env]] <- retainedInList
          }
        }
      }
      
      if(draw){
        if(!file.exists(paste("Data/intronRetention/", filename, "_nIn", nIn, "_Exp", cutoff, "_reThres", threshold, ".png", sep=""))){
          #DRAW PICTURE!!!!!
          qtl <- read.table(paste(location, filename, "_QTL.txt", sep=""), header=TRUE, row.names=1)
          
          png(file = paste("Data/intronRetention/", filename, "_nIn", nIn, "_Exp", cutoff, "_reThres", threshold, ".png", sep=""), bg="white", width=1024, height=1024)
          singlePlot(chr, rawexp, newexp, qtl, filename, intronID, exonExpID, menvironment, ind_tu, res, lodThreshold)
          dev.off()
          nPlot <- nPlot + 1
        } else{
          cat("Skipping", filename," because it exists\n")
        }
      } else{
        cat("Skipping", filename, "no retention\n")
      }
      
    }
  }
  
}


for(chr in 1:5){
  plotReBYchr(chr, menvironment, nIn=1, cutoff = 5, threshold = 7, lodThreshold = 4)
}
