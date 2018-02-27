#' Rscript for the snakemake pipeline for computing jump probabilities and including heatmaps
#' author Maryam Ghareghani

# install and load MaRyam package
library(devtools) 
install_git("git://github.com/friendsofstrandseq/MaRyam.git", branch = "master")
library("MaRyam")


hapProbsFile = snakemake@input[["haplotypeProbs"]]
GTprobsFile = snakemake@input[["genotypeProbs"]]

GTprobs <- read.table(GTprobsFile, stringsAsFactors = F, header = T)
hapProbs <- read.table(hapProbsFile, stringsAsFactors = F, header = T)

# getting the splitted prob tables including the jump probs
GTprobs.l.chrom <- addJumpProbs(GTprobs)
hapProbs.l.chrom <- addJumpProbs(hapProbs)

pdf(snakemake@output[[1]])
for (k in 1:length(GTprobs.l.chrom))
  lapply(1:length(GTprobs.l.chrom[[k]]), function(x) grid.arrange(plotHeatmapSegment(GTprobs.l.chrom[[k]][[x]])$heatmap.plt, plotHeatmapSegment(hapProbs.l.chrom[[k]][[x]])$heatmap.plt))
dev.off()
