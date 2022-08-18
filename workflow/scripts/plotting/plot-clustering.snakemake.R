log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = "message")
sink(file = log, type = "output")

source("workflow/scripts/plotting/plot-clustering.R")
plot.clustering(
    inputfile = snakemake@input[["sv_calls"]],
    bin.bed.filename = snakemake@input[["binbed"]],
    position.outputfile = snakemake@output[["position"]],
)
# chromosome.outputfile <- snakemake@output[["chromosome"]]

# args <- commandArgs(trailingOnly = TRUE)

# plot.clustering(
#     inputfile = args[1],
#     bin.bed.filename = args[2],
#     position.outputfile = args[3],
#     chromosome.outputfile = args[4]
# )
