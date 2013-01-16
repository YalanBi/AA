#
# Functions for analysing A. Thaliana Tiling Arrays
# last modified: 16-01-2013
# first written: 16-01-2013
# (c) 2013 GBIC Yalan Bi, Danny Arends, R.C. Jansen
#

#Idea---try to find the upper and lower separation points for intergenic
#key point1---arr[x,"gene"]:gene vs intergenic vs missing?
#							If arr[x,"gene"] == intergenic, intron vs intergenic?
#							If arr[x,"gene"] != intergenic, any intergenic between 2 genes?
#key point2---any intergenic left at the end?
makeSubsets <- function(filename, ){
  read.table(filename)########?????????#######
  previous_gene <- "" #characters
  previous_x <- 0 #numbers
  intergenic_seq <- NULL #a vector
  empty_seq <- NULL #a vector
  for(x in 1:nrow(arr)){
    #s <- proc.time()
    if(arr[x,"gene"] != "intergenic" && arr[x,"gene"] != ""){ #gene vs intergenic vs missing? ""---the missing one: without arr[x,"gene"] != "", if there are some intergenic between 2"", these intergenic can be regarded as introns, so will be missed in "inter_danny.txt".
      if(as.character(previous_gene) != as.character(arr[x,"gene"])){ #intron vs intergenic?
        cat(previous_x,as.character(previous_gene),"=",x, as.character(arr[x,"gene"]),"\n")
        if((x-1) == previous_x){ #any intergenic between 2 genes?
          #Too bad 2 genes without an intergenic region
        }else{
          #Found an REAL INTERGENIC probes
          for(y in seq((previous_x+1),(x-1),1)){ #with a for loop to check the missing part
            if(arr[y,"gene"] == "intergenic"){
              intergenic_seq <- c(intergenic_seq, y)
            }else{
              empty_seq <- c(empty_seq, y)
            }
          }
        }
      }
      previous_gene <- arr[x,"gene"]
      previous_x <- x
    }
    
    #cat(x,"took:",(proc.time()-s)[3],"\n")
  }
  x <- nrow(arr) #any intergenic left at the end?
  #AT THE END so see if we are having intergenic probes
  while(arr[x,"gene"] == "intergenic"){
    x <- (x-1)
  }
  #x is now at the line which has the last GENE !!!
  #to check whether while loop finds any intergenic new
  if(x != nrow(arr)){
    intergenic_seq <- c(intergenic_seq, seq(x+1, nrow(arr), 1))
  }

  write.table(arr[intergenic_seq,],file="inter_danny.txt")

  without <- arr[-intergenic_seq, ] #remove intergenic from the original data
  for(gene_name in unique(as.character(without[,"gene"]))){ #by unique(), each gene name check once
    if(gene_name != "intergenic" && gene_name != ""){
      ids <- which(as.character(without[,"gene"]) == gene_name)
      ids <- min(ids):max(ids) #contain the introns
      cat(gene_name,"\t",ids,"\n")
      write.table(without[ids, ], file=paste("gene_data1/", gene_name,".txt",sep=""))
    }else if(gene_name == ""){
      ids <- which(as.character(without[,"gene"]) == gene_name)
      cat("missing\t",ids,"\n")
      write.table(without[ids, ], file=paste("gene_data1/", "Missing.txt",sep=""))
    }
  }





  without <- arr[-intergenic_seq,] #remove intergenic from the original data
  for(gene_name in unique(as.character(without[,"gene"]))){ #by unique(), each gene name check once
    if(gene_name != "intergenic" && gene_name != ""){##########????????????##########
      ids <- which(as.character(without[,"gene"]) == gene_name)
      ids <- min(ids):max(ids) #contain the introns
      cat(gene_name,"\t",ids,"\n")
      write.table(without[ids, ], file=paste("gene_data/", gene_name,".txt",sep=""))
    }else if(gene_name == ""){########## ""never appears within a gene????????????##########
      ids <- which(as.character(without[,"gene"]) == gene_name)
      cat("missing\t",ids,"\n")
      write.table(without[ids, ], file=paste("gene_data/", "Missing.txt",sep=""))
    }
  }
}


