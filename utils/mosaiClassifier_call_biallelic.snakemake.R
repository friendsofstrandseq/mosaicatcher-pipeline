library(data.table)
source("utils/mosaiClassifier/makeSVcalls.R")

probs <- readRDS(snakemake@input[["probs"]])
llr <- as.numeric(snakemake@wildcards[["llr"]])

probs <- mosaiClassifierPostProcessing(probs)
probs <- forceBiallelic(probs)
tab <- makeSVCallSimple(probs, llr_thr = llr)

write.table(tab, file = snakemake@output[[1]], sep = "\t", quote = F, row.names = F, col.names = T)