log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

source('utils/haplotagProbs.R')

haplotagCounts <- fread(snakemake@input[["haplotag_table"]])
probs <- readRDS(snakemake@input[["sv_probs_table"]])

# FIXME: quick and dirty fix for off by one start coordinates of segments
haplotagCounts[,start:=start-1]

probs <- addHaploCountProbs(probs, haplotagCounts, alpha=0.05)

saveRDS(probs, file=snakemake@output[[1]])
