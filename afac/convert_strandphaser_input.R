# sink(snakemake@log[[1]])
library(data.table)
d = fread("/g/korbel2/weber/MosaiCatcher_output/Mosaicatcher_output_HGSVC_correct/segmentation/HTN7CAFXY_HG02492x02_19s003809-1-1/Selection_initial_strand_state")
print(d)
print(unique(d$cell))
print(unique(d$sample))
e = fread("/g/korbel2/weber/MosaiCatcher_output/Mosaicatcher_output_HGSVC_correct/counts/HTN7CAFXY_HG02492x02_19s003809-1-1/HTN7CAFXY_HG02492x02_19s003809-1-1.info")
tmp <- basename(e$bam)
e$bam <- strsplit(e$bam, split = ".", fixed = TRUE)[[1]][1]
print(e)
print(unique(e$cell))
print(unique(e$sample))
# f <- merge(d, e, by = c("sample", "cell"))[class == "WC", .(chrom, start, end, bam)]
# f <- merge(d, e, by = c("sample", "cell"))[class == "WC", .(chrom, start, end, bam)]
f <- merge(d, e, by = c("sample", "cell"))
print(f)
# write.table(f, file=snakemake@output[[1]], quote=F, row.names=F, col.names=F, sep="\t")