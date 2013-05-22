#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 18-04-2013
# first written: 18-04-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#*************************************************************** basic part ***************************************************************#
setwd("D:/Arabidopsis Arrays")


#*************************************************************** load part ***************************************************************#
#load exp genes
load(file="Data/fullModeMapping/expGenes.Rdata")


for(chr in 1:5){
  location <- paste0("Data/chr", chr, "_norm_hf_cor/")
  location_FM <- paste0("Data/fullModeMapping/chr", chr, "_norm_hf_cor/")
  
  genenames <- gsub(".txt", "", expGeneList[[chr]])
  chooseID <- sample(1:length(genenames), 100, replace = FALSE)
  cat("chr", chr, "copy genes are:", genenames[sort(chooseID)], "\n")
  
  for(filename in genenames[sort(chooseID)]){
    copyFiles <- dir(location)[grepl(filename, dir(location)) & !grepl("QTL", dir(location)) & !grepl(".png", dir(location))]
    file.copy(from = paste0(location, copyFiles), to = paste0("toDanny/rewExp/chr", chr, "_norm_hf_cor/", copyFiles))
    cat(filename, "copy", length(copyFiles), "files successfully!\n")
    
    copyFMfiles <- dir(location_FM)[grepl(filename, dir(location_FM))]
    file.copy(from = paste0(location_FM, copyFMfiles), to = paste0("toDanny/fullModeMapping/chr", chr, "_norm_hf_cor/", copyFMfiles))
    cat(filename, "copy", length(copyFMfiles), "FMfiles successfully!\n")
  }
}
