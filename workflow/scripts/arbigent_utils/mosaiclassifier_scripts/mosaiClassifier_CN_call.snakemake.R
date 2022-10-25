log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

#library(data.table)
source("workflow/scripts/arbigent_utils/mosaiclassifier_scripts/mosaiClassifier/makeSVcalls.R")

probs = readRDS(snakemake@input[["probs"]])

probs <- makeCNcall(probs)

write.table(probs, file = snakemake@output[[1]], sep = "\t", quote=F, row.names = F, col.names = T)
