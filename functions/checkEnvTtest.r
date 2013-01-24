#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 24-01-2013
# first written: 24-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

checkEnvTtest <- function(filename, location = "C:/Arabidopsis Arrays/Data/chr1/"){
  rawexp <- read.table(paste(location, filename, sep=""), header=TRUE, row.names=1)
  newexp <- rawexp[,16:163]
  s <- 1
  lgp <- NULL
  lgpmatrix <- NULL
  for(p in 1:nrow(rawexp)){
    env1 <- as.matrix(newexp[p,which(as.numeric(env[,2])==1)])
    env2 <- as.matrix(newexp[p,which(as.numeric(env[,2])==2)])
    env3 <- as.matrix(newexp[p,which(as.numeric(env[,2])==3)])
    env4 <- as.matrix(newexp[p,which(as.numeric(env[,2])==4)])
    env_mean <- mean(c(env1,env2,env3,env4))
    lgp <- -log10(c(t.test(env1, mu=env_mean)$p.value,t.test(env2, mu=env_mean)$p.value,t.test(env3, mu=env_mean)$p.value,t.test(env4, mu=env_mean)$p.value))
    #lgp <- c(lgp, round(lgp-0.5)+1) #####check for mycolor()!!!#####
    lgpmatrix <- rbind(lgpmatrix, lgp)
  }
  return(lgpmatrix)
}

checkEnvTtest("AT1G01010.txt")
