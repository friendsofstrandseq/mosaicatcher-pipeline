sink(snakemake@log[[1]])
library(data.table)
library(assertthat)
source("utils/mosaiClassifier/mosaiClassifier.R")


# Currently read files from the Snakemake pipeline
counts = fread(paste("zcat",snakemake@input[["counts"]]))
info   = fread(snakemake@input[["info"]])
strand = fread(snakemake@input[["states"]])
segs   = fread(snakemake@input[["bp"]])


# is there a normalization file given?
if ("norm" %in% names(snakemake@input) && length(snakemake@input[["norm"]])>0) {
  message("[MosaiClassifier] Read normalization from ", snakemake@input[["norm"]])
  normalization = fread(snakemake@input[["norm"]])
} else {
  normalization = NULL
}

# haplotypeMode?
if ("CW" %in% strand$class) {
  haplotypeMode = T
} else {
  haplotypeMode = F
}

d = mosaiClassifierPrepare(counts, info, strand, segs, normalization)
e = mosaiClassifierCalcProbs(d, maximumCN = 4, haplotypeMode = haplotypeMode)

saveRDS(e, file = snakemake@output[[1]])



