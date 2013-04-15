#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 12-04-2013
# first written: 12-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")


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

intronExp <- function(rawexp, intronID, threshould){
  for(p in intronID){
    for(i in 17:164){
      if(rawexp[p, i] >= threshould){
        rawexp[p, i] <- 1
      } else{
        rawexp[p, i] <- 0
      }
    }
  }
  return(rawexp[intronID, 17:164])
}


#*************************************************************** ratio, intron ***************************************************************#
threshould = 5

nGene <- NULL
nInall <- 0
nInHigh <- 0

for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  cat("\nNow chr", chr, "starts!\n")
  st <- proc.time()[3]
  counting <- 0
  
  #filename "AT1G01010.txt"
  for(filename in dir(location)[which(grepl(".txt", dir(location)) & !grepl("_QTL", dir(location)))]){
    rawexp <- read.table(paste0(location, filename), row.names=1, header=T)
    counting <- counting + 1
    
    probes_dir <- probesDir(rawexp)
    #cat(filename, "\nprobeDir:", probes_dir, "\n")
    #intronID <- introns of right direction
    intronID <- probes_dir[which(grepl("intron", rawexp[probes_dir, "tu"]))]
    #cat("introns:", intronID, "\n")
        
    if(length(intronID) > 0){
      #number of intron individuals higher than threshould
      nInHigh <- nInHigh + length(which(rawexp[intronID, 17:164] >= threshould))
      #number of all intron individuals
      nInall <- nInall + length(intronID)*148
    }
    
    #cat(counting, "\t")
  }
  et <- proc.time()[3]
  cat("chr", chr, "has", counting, "genes; and ends up in", et-st, "s\n")
  
  cat("until chr", chr, "there are", nInHigh, "intron individuals higher than", threshould, "and", nInall, "intron individuals in all\n")
  nGene <- c(nGene, counting)
}



hist(exonmean, breaks=50)
hist(intronmean, breaks=50, col=6, add=TRUE)





#************************************************************ intron vs exon ************************************************************#
setwd("D:/Arabidopsis Arrays")

intronGene <- NULL
for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  cat("Now chr", chr, "starts!\n")
  st <- proc.time()[3]
  
  genenames <- dir(location)[which(grepl(".txt", dir(location)) & !grepl("_QTL", dir(location)))]
  
  #filename "AT1G01010.txt"
  for(filename in genenames){
    rawexp <- read.table(paste0(location, filename), row.names=1, header=T)
    
    probes_dir <- probesDir(rawexp)
    
    #intronID <- introns of right direction
    intronID <- probes_dir[which(grepl("intron", rawexp[probes_dir, "tu"]))]
        
    if(length(intronID) > 0){
      intronGene <- c(intronGene, filename)
    }
  }
  et <- proc.time()[3]
  cat("intronGene finding on chr", chr, "finished in", et-st, "s\n")
}

set.seed(1)
selectGene <- intronGene[sample(1:length(intronGene), 1000, replace=F)]
set.seed(2)
selectGenee <- intronGene[sample(1:length(intronGene), 1000, replace=F)]
set.seed(3)
selectGeneee <- intronGene[sample(1:length(intronGene), 1000, replace=F)]

selectINmatrix <- NULL
selectExmatrix <- NULL
for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  for(fn in selectGene[which(grepl(paste0("AT", chr), selectGene))]){
    cat(fn, "\n")
    rawexp <- read.table(paste0(location, fn), row.names=1, header=T)
    
    probes_dir <- probesDir(rawexp)
    exonID <- probes_dir[which(grepl("tu", rawexp[probes_dir, "tu"]))]
    #intronID <- introns of right direction
    intronID <- probes_dir[which(grepl("intron", rawexp[probes_dir, "tu"]))]
    #cat("introns:", intronID, "\n")
    
    set.seed(1)
    selectINmatrix <- rbind(selectINmatrix, rawexp[intronID[sample(1:length(intronID), 1, replace=F)], 17:164])
    selectExmatrix <- rbind(selectExmatrix, rawexp[exonID[sample(1:length(exonID), 1, replace=F)], 17:164])
  }
  for(fn in selectGenee[which(grepl(paste0("AT", chr), selectGenee))]){
    cat(fn, "\n")
    rawexp <- read.table(paste0(location, fn), row.names=1, header=T)
    
    probes_dir <- probesDir(rawexp)
    exonID <- probes_dir[which(grepl("tu", rawexp[probes_dir, "tu"]))]
    #intronID <- introns of right direction
    
    set.seed(2)
    selectExmatrix <- rbind(selectExmatrix, rawexp[exonID[sample(1:length(exonID), 1, replace=F)], 17:164])
  }
  for(fn in selectGeneee[which(grepl(paste0("AT", chr), selectGeneee))]){
    cat(fn, "\n")
    rawexp <- read.table(paste0(location, fn), row.names=1, header=T)
    
    probes_dir <- probesDir(rawexp)
    exonID <- probes_dir[which(grepl("tu", rawexp[probes_dir, "tu"]))]
    #intronID <- introns of right direction
    
    set.seed(3)
    selectExmatrix <- rbind(selectExmatrix, rawexp[exonID[sample(1:length(exonID), 1, replace=F)], 17:164])
  }
}

hist(unlist(selectExmatrix), breaks=50)
hist(unlist(selectINmatrix), breaks=50, col=rgb(1,0,0.8,0.5), add=TRUE)




#************************************************************ intron vs intergenic ************************************************************#
setwd("D:/Arabidopsis Arrays")

#intronGene <- all gene containing intron probes of right direction
getIntronGene <- function(){
  intronGene <- NULL
  for(chr in 1:5){
    location <- paste0("Data/chr", chr, "_norm_hf_cor/")
    cat("Now chr", chr, "starts!\n")
    st <- proc.time()[3]
    
    genenames <- dir(location)[which(grepl(".txt", dir(location)) & !grepl("_QTL", dir(location)))]
    
    #filename "AT1G01010.txt"
    for(filename in genenames){
      rawexp <- read.table(paste0(location, filename), row.names=1, header=T)
      
      probes_dir <- probesDir(rawexp)
      
      #intronID <- introns of right direction
      intronID <- probes_dir[which(grepl("intron", rawexp[probes_dir, "tu"]))]
      
      if(length(intronID) > 0){
        intronGene <- c(intronGene, filename)
      }
    }
    et <- proc.time()[3]
    cat("intronGene finding on chr", chr, "finished in", et-st, "s\n")
  }
  return(intronGene)
}
intronGene <- getIntronGene()

#create random intron probes matrix
randomInmatrix <- function(intronGene, seed = 1){
  set.seed(seed)
  selectGene <- intronGene[sample(1:length(intronGene), 1000, replace=F)]

  selectINmatrix <- NULL
  for(chr in 1:5){
    location <- paste0("Data/chr", chr, "_norm_hf_cor/")
    for(fn in selectGene[which(grepl(paste0("AT", chr), selectGene))]){
      #cat(fn, "\n")
      rawexp <- read.table(paste0(location, fn), row.names=1, header=T)
      
      probes_dir <- probesDir(rawexp)
      exonID <- probes_dir[which(grepl("tu", rawexp[probes_dir, "tu"]))]
      #intronID <- introns of right direction
      intronID <- probes_dir[which(grepl("intron", rawexp[probes_dir, "tu"]))]
      #cat("introns:", intronID, "\n")
      
      set.seed(seed)
      selectINmatrix <- rbind(selectINmatrix, rawexp[intronID[sample(1:length(intronID), 1, replace=F)], 17:164])
    }
  }
  return(selectINmatrix)
}

#create random intergenic probes matrix
randomIntermatrix <- function(seed = 1){
  set.seed(seed)
  indIntergenic <- sample(1:96300, 200, replace=F)
  intergenicmatrix <- NULL

  for(chr in 1:5){
    cat("chr", chr, "starts!\n")
    filename <- paste0("Data/intergenic", chr, "_norm_hf_cor.txt")
    fp <- file(filename)
    open(fp)

    #cnt=0, this is the header line
    cnt <- 0
    mline <- readLines(fp, n=1)

    st <- proc.time()[3]
    while(length(mline) != 0L){
      cnt <- cnt + 1
      mline <- readLines(fp, n=1)
      if(cnt %in% indIntergenic){
        intergenicmatrix <- rbind(intergenicmatrix, as.numeric(strsplit(mline,"\t")[[1]][18:165]))
      }
    }
    close(fp)
    et <- proc.time()[3]
    cat("intergenic", chr, "finished in", et-st, "s!\n")
  }
  return(intergenicmatrix)
}

ratioIntron <- NULL
ratioInter <- NULL
for(seed in 1:5){
  cat(seed, "seed", "starts!\n")
  intronList <- unlist(randomInmatrix(intronGene, seed = seed))
  intergenicList <- unlist(randomIntermatrix(seed = seed))
  ratioIntron <- c(ratioIntron, length(which(intronList >= threshould))/length(intronList))
  ratioInter <- c(ratioInter, length(which(intergenicList >= threshould))/length(intergenicList))
  
  png(filename = paste0("intron-intergenic_seed", seed, ".png"), width = 1280, height = 960, bg = "white")
  hist(intronList, breaks=50, col=rgb(1,0,0.8,0.5))
  hist(intergenicList, breaks=50, col=rgb(0.05,0.75,0.35,0.35), add=T)
  dev.off()
}
ratioIntron
ratioInter
