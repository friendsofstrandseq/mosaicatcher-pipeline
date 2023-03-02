#log <- file(snakemake@log[[1]], open='wt')
#sink(file=log, type='message')
#sink(file=log, type='output')

source("utils/count_add_peakind.R")
count_add_peakind(inputfile = snakemake@input[["count_sort"]])

