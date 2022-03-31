sink(snakemake@log[[1]])
library(data.table)
qu = snakemake@params[["quantile"]]
print(qu)
d = fread(snakemake@input[[1]])
print(d)
# type = 1 is important to get discrete values!
d = d[, .SD[k == quantile(1:max(k), qu, type = 1)], by=chrom][,.(k, chrom, bps)]
print(d)
write.table(d, file = snakemake@output[[1]], row.names=F, quote=F, sep="\t")
