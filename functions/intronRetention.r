#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 22-03-2013
# first written: 22-03-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

setwd("D:/Arabidopsis Arrays")
menvironment <- read.table("Data/ann_env.txt",sep="\t")[,2]

ProbesSelect <- function(filname, exp_data){#check the expression of exons of right direction
  load(paste(location,"Classification_chr",chromosome,"_norm_hf_cor.Rdata",sep=""))
  class <- res[[gsub(".txt", "_QTL.txt", filename)]]#res[["AT1G01010_QTL.txt"]]
  if(unique(exp_data[,"strand"]) == "sense"){
    direction_id <- which(exp_data[, "direction"] == "reverse")
  }
  if(unique(exp_data[,"strand"]) == "complement"){
    direction_id <- which(exp_data[, "direction"] == "forward")
  }
  probesDirExp <- direction_id[which(direction_id %in% class$goodP)]
  return(probesDirExp)
}

#####Here!#####
checkIntronRetention <- function(chromosome = 1){
  location <- paste("D:/Arabidopsis Arrays/Data/chr", chromosome, "_norm_hf_cor/", sep="")
  for(filename in dir(location)[grepl(".txt", dir(location)) & !grepl("_QTL", dir(location))]){
    exp_data <- read.table(paste(location, filename, sep=""), header=T, row.names=1)
    probes_id <- ProbesSelect(filename, exp_data)
    for(env in 1:4){
      ind_env <- which(as.numeric(menvironment) == env)
    }
  }
}



tuORintron <- unique(aa[direction_id, "tu"])
for(a in 1:length(tuORintron)){
  tuORintron[a]
}