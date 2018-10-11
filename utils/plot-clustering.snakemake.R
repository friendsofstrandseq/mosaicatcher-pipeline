log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

source("utils/plot-clustering.R")
plot.clustering(inputfile = snakemake@input[["sv_calls"]], bin.bed.filename = snakemake@input[["binbed"]], position.outputfile = snakemake@output[["position"]], chromosome.outputfile = snakemake@output[["chromosome"]])
