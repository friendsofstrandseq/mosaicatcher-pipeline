log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = "message")
sink(file = log, type = "output")

system("LC_MEASUREMENT=C")

source("workflow/scripts/haplotagging_scripts/haplotagTable.R")

paired_end <- sub("\n", "", readChar(snakemake@input[["paired_end"]], file.info(snakemake@input[["paired_end"]])$size))

tab <- getHaplotagTable2(bedFile = snakemake@input[["bed"]], bam.file = snakemake@input[["bam"]], file.destination = snakemake@output[["tsv"]], paired_end = paired_end)

# SVplotting(snakemake@input[["sv_calls"]], snakemake@output[["barplot_high"]], snakemake@output[["barplot_med"]], snakemake@output[["barplot_low"]], snakemake@output[["barplot_rare"]])
