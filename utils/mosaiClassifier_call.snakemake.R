log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
source("utils/mosaiClassifier/makeSVcalls.R")

probs                = readRDS(snakemake@input[["probs"]])
llr                  = as.numeric(snakemake@wildcards[["llr"]])
use.pop.priors       = eval(parse(text=snakemake@wildcards[["pop_priors"]]))
regularizationFactor = 10^(-as.numeric(snakemake@wildcards[["regfactor"]]))

probs <- mosaiClassifierPostProcessing(probs, regularizationFactor = regularizationFactor)
tab <- makeSVCallSimple(probs, llr_thr = llr, use.pop.priors = use.pop.priors)

write.table(tab, file = snakemake@output[[1]], sep = "\t", quote=F, row.names = F, col.names = T)
