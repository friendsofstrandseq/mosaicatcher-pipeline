# log <- file(snakemake@log[[1]], open = "wt")
# sink(file = log, type = "message")
# sink(file = log, type = "output")
args <- commandArgs(trailingOnly = T)

source("utils/haplotagProbs.R")

haplotagCounts <- fread(args[1])
probs <- readRDS(args[2])

# FIXME : tmp solution to fix error : Segments must covered all bins, which happen for small scaffolds
chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

haplotagCounts <- haplotagCounts[haplotagCounts$chrom %in% chroms, ]


# FIXME: quick and dirty fix for off by one start coordinates of segments
haplotagCounts[, start := start - 1]

probs <- addHaploCountProbs(probs, haplotagCounts, alpha = 0.05)

saveRDS(probs, file = args[3])