library(data.table)
library(assertthat)
source("utils/mosaiClassifier/mosaiClassifier.R")


#counts = fread(paste("zcat","counts/simulation0-50000/50000_fixed.txt.gz"))
#info   = fread("counts/simulation0-50000/50000_fixed.info")
#strand = fread("strand_states/simulation0-50000/final.txt")
#segs   = fread("segmentation2/simulation0-50000/50000_fixed.medium.txt")


# Currently read files from the Snakemake pipeline
counts = fread(paste("zcat",snakemake@input[["counts"]]))
info   = fread(snakemake@input[["info"]])
strand = fread(snakemake@input[["states"]])
segs   = fread(snakemake@input[["bp"]])


d = mosaiClassifierPrepare(counts, info, strand, segs)
d = mosaiClassifierCalcProbs(d)
d = mosaiClassifierPostProcessing(d)
saveRDS(d, file = snakemake@output[[1]])



