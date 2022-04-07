log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = "message")
sink(file = log, type = "output")

source("utils/haplotagTable.R")


tab <- getHaplotagTable2(bedFile = snakemake@input[["bed"]], bam.file = snakemake@input[["bam"]], file.destination = snakemake@output[["tsv"]])

# SVplotting(snakemake@input[["sv_calls"]], snakemake@output[["barplot_high"]], snakemake@output[["barplot_med"]], snakemake@output[["barplot_low"]], snakemake@output[["barplot_rare"]])