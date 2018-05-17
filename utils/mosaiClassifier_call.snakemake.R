library(data.table)
source("utils/mosaiClassifier/makeSVcalls.R")

probs = readRDS(snakemake@input[["probs"]])
llr   = as.numeric(snakemake@params[["llr"]])

tab <- makeSVCallSimple(probs)

write.table(tab, file = snakemake@output[[1]], sep = "\t", quote=F, row.names = F, col.names = T)
