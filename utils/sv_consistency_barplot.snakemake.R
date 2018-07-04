log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

source("utils/sv_consistency_barplot.R")

SVplotting(snakemake@input[["sv_calls"]], snakemake@output[["barplot_high"]], snakemake@output[["barplot_med"]], snakemake@output[["barplot_low"]], snakemake@output[["barplot_rare"]])
