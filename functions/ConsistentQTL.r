

setwd("D:/Arabidopsis Arrays")





consistentQTL <- function(){
  
}


for(chr in 1:4){
  location <- paste0("Data/newTU/geneObject_chr", chr, "/")
  for(name in dir(location)[grepl(".R", dir(location))]){
    load(paste(location, name, sep=""))
    
  }
  
}
