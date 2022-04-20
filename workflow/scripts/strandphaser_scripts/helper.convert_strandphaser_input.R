sink(snakemake@log[[1]])
library(data.table)

d <- fread(snakemake@input[["states"]])

e <- fread(snakemake@input[["info"]])
e$bam <- basename(e$bam)
f <- merge(d, e, by = c("sample", "cell"))[class == "WC", .(chrom, start, end, bam)]

write.table(f, file = snakemake@output[[1]], quote = F, row.names = F, col.names = F, sep = "\t")