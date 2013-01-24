#
#
# Some code: set PATH=%PATH%;C:\Documents and Settings\S2347261\Git\cmd
#
#

plotProbes <- function(filename, location = "C:/Arabidopsis Arrays/Data/chr1/"){
  qtl <- read.table(paste(location, filename, sep=""), header=TRUE, row.names=1)
  mv <- max(qtl)
  plot(c(0,ncol(qtl)), c(0,mv), t='n')
  for(x in 1:nrow(qtl)){
    points(as.numeric(qtl[x,]), t='l', col=x)
  }
}

plotProbes("AT1G01470_QTL.txt")