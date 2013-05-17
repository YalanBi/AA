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

#5'/3' AS test
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
    
    #for 5'/3' AS, at least 2 exons in a gene!!!
    if(length(uniqueExon) >= 2){
      #cat("we have >= 2 exons!\n")
      
      for(testExon in uniqueExon[c(1, length(uniqueExon))]){
        #ind <- judge which probe in exonID is of 1st exon name (T/F)
        ind <- rawexp[exonID, "tu"] == testExon
        #cat(as.character(testExon), "has probes", exonID[ind], "\n")
        
        #at least 6 probes in a group, try to avoid this case---one is highly expressed, the other is lowly expressed...
        if(length(which(ind)) >= 6){
          #cat("I'm", as.character(testExon), ">= 6 good probes, t.test me and remember my gene name!\n")
          
          #to find max difference between each probe in exp, to group them
          dff <- NULL
          for(p in exonID[ind][2:length(exonID[ind])]){
            dff <- c(dff, abs())
          }
        }
      }
    #*************************first, test the 1st exon/5' exon*************************#
      indF <- rawexp[exonID, "tu"] == uniqueExon[1]
      
      
      
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
