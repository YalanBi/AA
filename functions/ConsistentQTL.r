

setwd("D:/Arabidopsis Arrays")




consistentQTL <- function(){
  
}


#load QTL files through the entire file
loadQTLfile <- function(chr = 1){
  x <- 0
  matrix <- NULL
  mnames <- NULL
  fp <- file(paste0("Data/newTU/genesByTU_chr", chr, "_norm_hf_cor_QTL_t.txt"))
  open(fp)
  line <- readLines(fp, 1)
  tnames <- strsplit(line,"\t")[[1]][-1]
  line <- readLines(fp, 1)
  x <- x+1
  
  while(length(line) > 0 && line != ""){
    data <- strsplit(line,"\t")[[1]]
    mnames <- c(mnames, data[1])
    cat("Done", x, "\n")
    
    data <- round(as.numeric(data[-1]), d=2)
    matrix <- cbind(matrix, data)
    line <- readLines(fp, 1)
    
    x <- x+1
  }
  colnames(matrix) <- mnames
  rownames(matrix) <- tnames
  
  close(fp)
  invisible(matrix)
}

decideLoading <- function(chr){
  if(chr == 1 | chr == 5){
    qtlchr <- loadQTLfile(chr = chr)
  }else{
    qtlchr <- read.table(paste("Data/newTU/genesByTU_chr", chr, "_norm_hf_cor_QTL.txt", sep=""), row.names=1, header=TRUE)
  }
}

qtlchr <- decideLoading(chr)


#load QTL files through splitted files
for(chr in 1:5){
  location <- paste0("Data/newTU/geneObject_chr", chr, "/")
  qtlchr <- NULL
  for(name in dir(location)[grepl(".R", dir(location))]){
    load(paste(location, name, sep=""))
    qtlmatrix <- NULL
    for(tu in 2:length(res)){
      qtl <- res[tu][[1]]$QTL
      qtlmatrix <- rbind(qtlmatrix, qtl)
      
      
    }
    rownames(qtlmatrix) <- rownames(res$Exp)
    qtlchr <- rbind(qtlchr, qtlmatrix)
  }

}







peak.detect <- function(qtlprofiles, verbose = FALSE){
  cutoff = -log10(0.05 / (716 * nrow(qtlprofiles)))
  if(verbose) cat("Starting peak detection LOD >=", cutoff, "\n")
  mmatrix <- NULL
  for(x in 1:nrow(qtlprofiles)){
    peak <- FALSE
    curmax <- 0
    curmaxindex <- 1
    marker <- 1
    maximums <- NULL
    mrow <- rep(0,ncol(qtlprofiles))
    for(ab in (qtlprofiles[x,]>cutoff | qtlprofiles[x,]<(-cutoff))){
      if(ab){
        peak <- TRUE
        if(qtlprofiles[x,marker]/abs(qtlprofiles[x,marker]) > 0){
          if(qtlprofiles[x,marker] > curmax){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- marker
          }
        }else{
          if(qtlprofiles[x,marker] < (-curmax)){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- -marker
          }
        }
        if(ncol(qtlprofiles)==marker){
          if(curmax!=0) maximums <- c(maximums,curmaxindex)
        }
      }else{
        if(curmax!=0) maximums <- c(maximums,curmaxindex)
        peak <- FALSE
        curmax <- 0
      }
      marker <- marker+1
    }
    mrow[which(qtlprofiles[x,] > cutoff)] <- 1
    mrow[which(qtlprofiles[x,] < -cutoff)] <- -1
    for(a in which(maximums>0)){
      mrow[maximums[a]] <- 2
    }
    for(b in which(maximums<0)){
      mrow[(-maximums[b])] <- -2
    }
    mmatrix <- rbind(mmatrix,mrow)
  }
  rownames(mmatrix) <- rownames(qtlprofiles)
  colnames(mmatrix) <- colnames(qtlprofiles)
  #mmatrix
  write.table(mmatrix, file="Data/newTU/genesByTU_chr1_norm_hf_cor_peak.txt")#change output filename
}

aa <- peak.detect(qtlprofiles = qtlchr, verbose = TRUE)
