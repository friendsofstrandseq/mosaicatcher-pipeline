library(ggplot2)
source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")
library('ggtree')
library("ape")

dir <- "/home/maryam/research/hackathons/C7-7Jun"

tab <- data.table::fread(file.path(dir, "BFB_cell_CNs.table"))
CN.mat <- split(tab$CN, tab$cell) 
CN.mat <- do.call(rbind, CN.mat)
euc.dist <- dist(CN.mat, method = "euclidean")
hc <- hclust(euc.dist, method = "ward.D2")
plot(hc)

plot(hclust(CN.mat, ))
ggtree(as.phylo(hc)) + facet_grid()

ggplot(data=cellCNcalls)