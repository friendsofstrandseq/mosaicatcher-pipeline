log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = "message")
sink(file = log, type = "output")

source("workflow/scripts/plotting/sv_consistency_barplot.R")

SVplotting(inputfile = snakemake@input[["sv_calls"]], outputfile.byPOS = snakemake@output[["barplot_bypos"]], outputfile.byVAF = snakemake@output[["barplot_byaf"]])