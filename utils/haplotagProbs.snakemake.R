sink(snakemake@log[[1]])

haplotagCounts <- snakemake@input[[haplotag_table]]
probs <- snakemake@input[[sv_probs_table]]

probs <- addHaploCountProbs(probs, haploCounts, alpha)

save(probs, file=snakemake@output[[1]])