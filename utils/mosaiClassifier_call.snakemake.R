log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
source("utils/mosaiClassifier/makeSVcalls.R")

probs                = readRDS(snakemake@input[["probs"]])
llr                  = as.numeric(snakemake@wildcards[["llr"]])
use.pop.priors       = eval(parse(text=snakemake@wildcards[["pop_priors"]]))
use.haplotags        = eval(parse(text=snakemake@wildcards[["use_haplotags"]]))
regularizationFactor = 10^(-as.numeric(snakemake@wildcards[["regfactor"]]))
genotype.cutoff      = as.numeric(snakemake@wildcards[["gtcutoff"]])
minFrac.used.bins    = as.numeric(snakemake@params[["minFrac_used_bins"]])
bin.size     	     = as.numeric(snakemake@wildcards[["window"]])

probs <- mosaiClassifierPostProcessing(probs, regularizationFactor = regularizationFactor)
tab <- makeSVCallSimple(probs, llr_thr = llr, use.pop.priors = use.pop.priors, use.haplotags = use.haplotags, genotype.cutoff = genotype.cutoff, bin.size, minFrac.used.bins = minFrac.used.bins)

write.table(tab, file = snakemake@output[[1]], sep = "\t", quote=F, row.names = F, col.names = T)
