library(data.table)
source("/g/korbel2/weber/workspace/mosaicatcher-update/workflow/scripts/mosaiclassifier_scripts/mosaiClassifier/makeSVcalls.R")

probs <- readRDS("/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG_dryrun/mosaiclassifier/haplotag_likelihoods/RPE-BM510.Rdata")
llr <- as.numeric(4)
use.pop.priors <- eval(parse(text = TRUE))
use.haplotags <- eval(parse(text = TRUE))
regularizationFactor <- 10^(-as.numeric(6))
genotype.cutoff <- as.numeric(0)
minFrac.used.bins <- as.numeric(0.8)
bin.size <- as.numeric(100000)

print(probs)

print(llr)

probs <- mosaiClassifierPostProcessing(probs, regularizationFactor = regularizationFactor)

print(probs)

tab <- makeSVCallSimple(probs, llr_thr = llr, use.pop.priors = use.pop.priors, use.haplotags = use.haplotags, genotype.cutoff = genotype.cutoff, bin.size, minFrac.used.bins = minFrac.used.bins)

print(tab)


# write.table(tab, file = , sep = "\t", quote = F, row.names = F, col.names = T)