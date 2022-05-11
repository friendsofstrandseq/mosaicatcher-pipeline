library(data.table)
qu = 0.4
print(qu)
d = fread("/g/korbel2/weber/MosaiCatcher_output/segmentation/ERR2940607/100000.txt")
print(d)
print('\n')
print(1:max(70))
print(quantile(1:max(70), qu, type = 1))
stop()
# type = 1 is important to get discrete values!
d = d[, .SD[k == quantile(1:max(k), qu, type = 1)], by=chrom][,.(k, chrom, bps)]
print(d)
# write.table(d, file = snakemake@output[[1]], row.names=F, quote=F, sep="\t")