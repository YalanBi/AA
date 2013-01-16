#P6-New map check

setwd("X:/Data/Desktop/Arabidopsis Arrays")
newgenotype <- read.table("refined map/genotypes.txt")
dim(newgenotype)
#[1] 148 716
newmap <- read.table("refined map/map.txt")
dim(newmap)
#[1] 716   2
genotypes_n <- read.table("Data/genotypes_n.txt")
dim(genotypes_n)
#[1] 164  69
rownames(genotypes_n)[which(! rownames(genotypes_n) %in% rownames(newgenotype))]
#[1] "RIL59" "RIL60" "RIL61" "RIL62" "RIL63" "RIL64" "RIL65" "RIL66" "RIL67"
#[10] "RIL68" "RIL69" "RIL70" "RIL71" "RIL72" "RIL73" "RIL74"
