sink(snakemake@log[[1]])
library(data.table)
source("utils/mosaiClassifier/makeSVcalls.R")

probs = readRDS(snakemake@input[["probs"]])
llr   = as.numeric(snakemake@params[["llr"]])

probs <- mosaiClassifierPostProcessing(probs)
probs <- forceBiallelic(probs)
tab <- makeSVCallSimple(probs)

write.table(tab, file = snakemake@output[[1]], sep = "\t", quote=F, row.names = F, col.names = T)
