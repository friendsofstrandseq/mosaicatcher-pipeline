log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = "message")
sink(file = log, type = "output")

library(data.table)
source("workflow/scripts/mosaiclassifier_scripts/mosaiClassifier/makeSVcalls.R")

# probs <- readRDS(snakemake@input[["probs"]])
# llr <- as.numeric(snakemake@wildcards[["llr"]])
# use.pop.priors <- eval(parse(text = snakemake@wildcards[["pop_priors"]]))
# use.haplotags <- eval(parse(text = snakemake@wildcards[["use_haplotags"]]))
# regularizationFactor <- 10^(-as.numeric(snakemake@wildcards[["regfactor"]]))
# genotype.cutoff <- as.numeric(snakemake@wildcards[["gtcutoff"]])
# minFrac.used.bins <- as.numeric(snakemake@params[["minFrac_used_bins"]])
# bin.size <- as.numeric(snakemake@params[["window"]])


probs <- readRDS(snakemake@input[["probs"]])
llr <- as.numeric(snakemake@params[["llr"]])
use.pop.priors <- eval(parse(text = snakemake@params[["pop_priors"]]))
use.haplotags <- eval(parse(text = snakemake@params[["use_haplotags"]]))
regularizationFactor <- 10^(-as.numeric(snakemake@params[["regfactor"]]))
genotype.cutoff <- as.numeric(snakemake@params[["gtcutoff"]])
minFrac.used.bins <- as.numeric(snakemake@params[["minFrac_used_bins"]])
bin.size <- as.numeric(snakemake@params[["window"]])

# print(probs)
# print(llr)
# print(use.pop.priors)
# print(regularizationFactor)
# print(genotype.cutoff)
# print(minFrac.used.bins)
# print(bin.size)
# stop()

probs <- mosaiClassifierPostProcessing(probs, regularizationFactor = regularizationFactor)

# print(probs)

tab <- makeSVCallSimple(probs, llr_thr = llr, use.pop.priors = use.pop.priors, use.haplotags = use.haplotags, genotype.cutoff = genotype.cutoff, bin.size, minFrac.used.bins = minFrac.used.bins)

# print(tab)

write.table(tab, file = snakemake@output[[1]], sep = "\t", quote = F, row.names = F, col.names = T)