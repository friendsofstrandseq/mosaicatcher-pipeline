#' Rscript for the snakemake pipeline for computing jump probabilities and including heatmaps
#' author Maryam Ghareghani

sink(snakemake@log[[1]])

.libPaths( c( snakemake@params[["r_package_path"]],.libPaths()) )
suppressPackageStartupMessages(library(MaRyam))

hapProbsFile = snakemake@input[["haplotypeProbs"]]
GTprobsFile = snakemake@input[["genotypeProbs"]]

GTprobs <- read.table(GTprobsFile, stringsAsFactors = F, header = T)
hapProbs <- read.table(hapProbsFile, stringsAsFactors = F, header = T)

# Converting the prob tables to GT classe
GTclassProbs <- get_GT_class_prob_table(GTprobs)
hapClassProbs <- get_GT_class_prob_table(GTprobs)

# getting the splitted prob tables including the jump probs
GTprobs.l.chrom <- addJumpProbs(GTclassProbs)
hapProbs.l.chrom <- addJumpProbs(hapClassProbs)

pdf(snakemake@output[[1]])
for (k in 1:length(GTprobs.l.chrom))
  lapply(1:length(GTprobs.l.chrom[[k]]), function(x) grid.arrange(plotHeatmapSegment(GTprobs.l.chrom[[k]][[x]]), plotHeatmapSegment(hapProbs.l.chrom[[k]][[x]])))
dev.off()
