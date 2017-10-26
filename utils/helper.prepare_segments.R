library(data.table)
qu = snakemake@params[["quantile"]]
d = fread(snakemake@input[[1]])
# type = 1 is important to get discrete values!
d = d[, .SD[k == quantile(1:max(k), qu, type = 1)], by=chrom][,.(k, chrom, bps = breakpoint)]
write.table(d, file = snakemake@output[[1]], row.names=F, quote=F, sep="\t")