sink(snakemake@log[[1]])
library(data.table)
library(assertthat)
source("scripts/mosaiclassifier_scripts/mosaiClassifier/mosaiClassifier.R")


# Currently read files from the Snakemake pipeline
counts <- fread(paste("zcat", snakemake@input[["counts"]]))
info <- fread(snakemake@input[["info"]])
strand <- fread(snakemake@input[["states"]])
segs <- fread(snakemake@input[["bp"]])


# FIXME : tmp solution to fix error : Segments must covered all bins, which happen for small scaffolds
chroms <- snakemake@config[["chromosomes"]]

counts <- counts[counts$chrom %in% chroms, ]
strand <- strand[strand$chrom %in% chroms, ]
segs <- segs[segs$chrom %in% chroms, ]

# haplotypeMode?
if ("CW" %in% strand$class) {
  haplotypeMode <- T
} else {
  haplotypeMode <- F
}

d <- mosaiClassifierPrepare(counts, info, strand, segs)
e <- mosaiClassifierCalcProbs(d, maximumCN = 4, haplotypeMode = haplotypeMode)

saveRDS(e, file = snakemake@output[[1]])