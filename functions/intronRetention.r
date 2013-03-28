#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 26-03-2013
# first written: 22-03-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#
###########
# library #
###########

setwd("D:/Arabidopsis Arrays")
menvironment <- read.table("Data/ann_env.txt", sep="\t")[,2]

#############
# Functions #
#############

#####This function does more stuff#####
ProbesSelect <- function(location, chromosome, filename, exp_data){#check the expression of exons of right direction
  load(paste(location, "Classification_chr", chromosome, "_norm_hf_cor.Rdata", sep=""))
  class <- res[[gsub(".txt", "_QTL.txt", filename)]]#res[["AT1G01010_QTL.txt"]]
  probes <- NULL
  if(unique(exp_data[,"strand"]) == "sense"){
    direction_id <- which(exp_data[, "direction"] == "reverse")
  }
  if(unique(exp_data[,"strand"]) == "complement"){
    direction_id <- which(exp_data[, "direction"] == "forward")
  }
  probesDirExp <- direction_id[which(direction_id %in% class$goodP)]
  probes$In <- probesDirExp[which(probesDirExp %in% class$introP)]
  probes$Ex <- probesDirExp[-which(probesDirExp %in% class$introP)]
  return(probes)
}

isExpressed <- function(LODscore){#-log10(0.05/100000)=6.30
  if(LODscore >= 7) return(0) #0 means significant difference, not exp
  return(1) #1 means unsignificant difference, exp!!!
}


#####This function..... does stuff#####
checkIntronRetention <- function(chromosome = 1, startAt = "AT1G01010"){
  location <- paste("D:/Arabidopsis Arrays/Data/chr", chromosome, "_norm_hf_cor/", sep="")
  intronmatrix <- NULL
  gene_intron <- NULL
  allfiles <- dir(location)[grepl(".txt", dir(location)) & !grepl("_QTL", dir(location))]
  s <- which(grepl(startAt, allfiles))
  e <- length(allfiles)
  for(filename in allfiles[s:e]){
    cat(filename, "loading...\n")
    exp_data <- read.table(paste(location, filename, sep=""), header=T, row.names=1)
    new_exp <- exp_data[,17:ncol(exp_data)]
    probes_id <- ProbesSelect(location, chromosome, filename, exp_data)
    
    #For every intron, do stuff 
    for(i in unique(exp_data[,"tu"][probes_id$In])){
      intronline <- NULL
      intron_id <- probes_id$In[which(probes_id$In %in% which(exp_data[,"tu"] == i))]
      #For every environment you do ttest
      for(env in 1:4){
        ind_env <- which(as.numeric(menvironment) == env)
        #What is in Y?
        y <- as.numeric(unlist(new_exp[probes_id$Ex, ind_env]))
        #What is in X?
        x <- as.numeric(unlist(new_exp[intron_id, ind_env]))
        #T test on what?
        if(length(x) > 2 && length(y) > 2){
          tres <- t.test(x, y)
          LODscore <- -log10(tres$p.value)
          expIn <- isExpressed(LODscore)
          intronline <- c(intronline, round(tres$estimate[1], d=2), round(tres$statistic,d=3), round(LODscore,d=1), expIn)
        }
      }
      if(!is.null(intronline)){
        gene_intron <- c(gene_intron, paste(gsub(".txt", "", filename), "_", i, sep=""))
        intronmatrix <- rbind(intronmatrix, intronline)
      }
    }
  }
  rownames(intronmatrix) <- gene_intron
  colnames(intronmatrix) <- c("Env1_m", "Env1_t", "Env1_p", "Env1_Exp", "Env2_m", "Env2_t", "Env2_p", "Env2_Exp", "Env3_m", "Env3_t", "Env3_p", "Env3_Exp", "Env4_m", "Env4_t", "Env4_p", "Env4_Exp")
  write.table(intronmatrix, file=paste("Data/intronRetention_chr", chromosome, ".txt", sep=""))
  invisible(intronmatrix)
}

########
# Code #
########
checkIntronRetention(chromosome = 1, startAt = "AT1G01010")
