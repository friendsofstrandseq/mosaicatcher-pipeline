# sink(snakemake@log[[1]])
library(data.table);
d = fread("/g/korbel2/weber/MosaiCatcher_output_sample_HJ/strand_states/TALL03-DEA5/100000.selected_j0.1_s0.5_scedist20/initial_strand_state")
print(d)
e = fread("/g/korbel2/weber/MosaiCatcher_output_sample_HJ/counts/TALL03-DEA5/100000.info")
tmp = basename(e$bam)
e$bam = strsplit(e$bam, split=".", fixed=TRUE)[[1]][1]
print(e)
f = merge(d, e, by = c("sample","cell"))[class == "WC", .(chrom,start,end,bam)]
# f = merge(d, e, by = c("sample","cell"))
print(f)
# write.table(f, file=snakemake@output[[1]], quote=F, row.names=F, col.names=F, sep="\t")
