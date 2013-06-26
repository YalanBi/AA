#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 15-04-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#


#*********************************************** this is the final version for spliting original file into small pieces ^_^ **********************************************#
setwd("D:/Arabidopsis Arrays")

splitChromosome <- function(filename, chr=1){
  st <- proc.time()
  
  fp <- file(filename)
  open(fp)
  cnt <- 0
  mline <- ""
  previous_gene <- ""
  previous_x <- 0
  intergenic_seq <- NULL
  empty_seq <- NULL
  header <- c("ID", strsplit(readLines(fp, n=1), "\t")[[1]])
  igfilename <- paste0("Data/Raw/intergenic", chr, "_norm_hf_cor.txt")
  cat(paste(header, collapse="\t"), "\n", file=igfilename, sep="")
  mline <- readLines(fp, n=1)
  
  while(length(mline) != 0L){
    elements <- strsplit(mline, "\t")[[1]]
    names(elements) <- header
    
    # Everything done we can start analysis of this probe !
    if(elements["gene"] != "intergenic" && elements["gene"] != ""){ # A gene !
      if(as.character(previous_gene) != as.character(elements["gene"])){ # A new gene
        cat(previous_x, as.character(previous_gene), "=", cnt, as.character(elements["gene"]), "in", (et-st)[3], "secs\n")
        outfile <- paste0("Data/Raw/chr", chr, "_norm_hf_cor/", as.character(elements["gene"]), ".txt")
        cat(paste(header,collapse="\t"), "\n", file=outfile)
      } #Write the probe to the gene file
      cat(mline, "\n", file=outfile, sep="", append=TRUE)
      previous_gene <- elements["gene"]
      previous_x <- cnt
    }else{ # Write the probe to the intergenic file
      cat(mline, "\n", file=igfilename, sep="", append=TRUE)
    }
    cnt <- cnt + 1
    et <- proc.time()
    mline <- readLines(fp, n=1)
    #cat(cnt, "element took: ", (et-st)[3]," seconds\n")
  }
  close(fp)
  cat(filename, "contains", cnt, "probes\n")
}
#splitChromosome("Yalan map/offspring_phenotypes_chr2_norm_hf_cor.txt", chr = 2)


for(filename in dir("Yalan map/")[grepl("cor",dir("Yalan map/"))]){
  chr <- which(filename == dir("Yalan map/")[grepl("cor",dir("Yalan map/"))])
  splitChromosome(paste("Yalan map/", filename, sep=""), chr)
}
