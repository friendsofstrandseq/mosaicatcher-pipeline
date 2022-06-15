library(data.table)
source("workflow/scripts/mosaiclassifier_scripts/mosaiClassifier/makeSVcalls.R")

probs <- readRDS(snakemake@input[["probs"]])
llr <- as.numeric(snakemake@wildcards[["llr"]])
bin_size <- as.numeric(snakemake@wildcards[["window"]])

probs <- mosaiClassifierPostProcessing(probs)
probs <- forceBiallelic(probs)
tab <- makeSVCallSimple(probs, llr_thr = llr, bin.size = bin_size)

write.table(tab, file = snakemake@output[[1]], sep = "\t", quote = F, row.names = F, col.names = T)