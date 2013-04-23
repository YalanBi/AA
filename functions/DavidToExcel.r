#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 15-04-2013
# first written: 15-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")


#************************************************************* split&save part *************************************************************#
raw <- read.table("clusterBlack_David.txt",fill=T,sep="\t")
selected <- raw[,c(1,2,3,5,6,12)]
n <- length(grep("Annotation Cluster",selected[,1]))
breaks <- c(grep("Annotation Cluster",selected[,1]), nrow(selected)+1)
data <- vector("list",n)
for(x in 1:n){
  data[[x]] <- selected[breaks[x]:(breaks[x+1]-1),]
}

for(x in 1:length(data)){
  data[[x]] <- data[[x]][-1, ]
  colnames(data[[x]]) <- as.character(unlist(data[[x]][1, ]))
  data[[x]] <- data[[x]][-1, ]
  data[[x]][,3] <- as.numeric(as.character(data[[x]][,3]))
  write.table(data[[x]], file=paste("clusterBlack_cluster", x, ".txt", sep=""), quote=FALSE, sep='\t')
}
