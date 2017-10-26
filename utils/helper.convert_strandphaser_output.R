library(data.table);
library(assertthat)
d = fread(snakemake@input[["phased_states"]])
e = fread(snakemake@input[["info"]])
g = fread(snakemake@intput[["initial_states"]])

d$bam = basename(d$bam);

e$bam = e$cell
e$cell = NULL
e$sample = NULL
f = merge(d, e, by = "bam")[, .(chrom,start,end,sample,cell,class)]

# Note that there is still a bug in Venla's strand state detection.
g = merge(g,f, by = c("chrom","start","end","sample","cell"), all.x = T)

# Overwrite with David's phased strand state if available!
g = g[, class := ifelse(!is.na(class.y), class.y, class.x)][]
g$class.x = NULL
g$class.y = NULL
g = g[, .(chrom, start, end, sample, cell, class)]

write.table(g, file=snakemake@output[[1]], quote=F, row.names=F, col.names=F, sep="\t")