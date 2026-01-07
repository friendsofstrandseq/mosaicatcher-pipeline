library(ggplot2)
library(reshape2)

# cellwise <- "cell" %in% names(snakemake@wildcards)
if (("cell" %in% names(snakemake@wildcards)) == TRUE) {
    wc_cell_row_plate <- snakemake@wildcards[["cell"]]
    type <- "Cell"
} else if(("row" %in% names(snakemake@wildcards)) == TRUE) {
    wc_cell_row_plate <- snakemake@wildcards[["row"]]
    type <- "Row"
} else {
    wc_cell_row_plate <- "Full plate"
    type <- "Type"
}

# zcat out.tsv.gz  | grep "^GC" > scripts/gc.table

args = commandArgs(trailingOnly=TRUE)
# x = read.table(args[1], header=T)
x = read.table(snakemake@input[["table"]], header = T)
df = dcast(x, GCcontent ~ Sample)
refcol = (1:ncol(df))[colnames(df) == "Reference"]
if (refcol == 2){gc_mean = sum(df[,1]*df[,3])}
if (refcol == 3){gc_mean = sum(df[,1]*df[,2])}

if (refcol == 3) { datacol = 2; } else { datacol = 3; }
df$observed = df[,datacol] / df[,refcol]
df$expected = 1.0
df = df[df$GCcontent > 0.2 & df$GCcontent < 0.7,]
sse = sum((df$observed - df$expected)^2)
df = melt(df, id.vars=colnames(df[,1:3]))

# png(args[2])
png(snakemake@output[["gcdist_plot"]])
p = ggplot(data=x, aes(x=GCcontent, y=fractionOfReads))
p = p + geom_freqpoly(aes(color=Sample), stat="identity")
p = p + xlab("GC content")
p = p + ylab("Normalized count")
p = p + ggtitle(paste0("Sample = ", snakemake@wildcards[["sample"]], "\n", type, "= ", wc_cell_row_plate, "\nGC content mean = ", format(round(gc_mean, 2), nsmall = 2))) + theme(plot.title = element_text(size = 12))

p
dev.off()

# png(args[3])
png(snakemake@output[["gcdevi_plot"]])
p = ggplot(data=df, aes(x=GCcontent, y=value)) + geom_line(aes(color=variable))
p = p + xlab("GC content")
p = p + ylab("Sample to reference")
p = p + ggtitle(paste0("Sample = ", snakemake@wildcards[["sample"]], "\n", type, "= ", wc_cell_row_plate, "\nSSE = ", format(round(sse, 2), nsmall = 2))) + theme(plot.title = element_text(size = 12))
p
dev.off()