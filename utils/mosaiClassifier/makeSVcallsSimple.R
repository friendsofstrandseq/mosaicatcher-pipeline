sink(snakemake@log[[1]])
library(data.table)
library(assertthat)

f_maryam = snakemake@input[["probs"]]
f_maryam
f_info   = snakemake@input[["info"]]
f_info
sample   = snakemake@params[["sample_name"]]
sample
f_bamNames = snakemake@input[["bamNames"]]
f_bamNames

d = fread(f_maryam, header = T)
print("Maryam's data:")
d
f = fread(f_info)
print("Info data:")
f
b = fread(f_bamNames)
b = b[order(cell_id),]
print("bamNames data:")
b

# Map cell number to cell name
assert_that(max(d$cells) <= nrow(f))
assert_that(max(d$cells) <= max(b$cell_id))
assert_that(all(1:nrow(b) == b$cell_id))    # Make sure that bamNames are sorted and exactly match numbers 1...n 

d$cell = b[d$cells,]$cell_name
d$sample = sample
d

# Get log likelihood ratio of major SV classes
d <- d[, .(chrom = chr,
           start = format(start, scientific=F),
           end   = format(end, scientific=F),
           sample, 
           cell, 
           type = toupper(types),
           p_del_hom = log(`0000`) - log(`1010`), 
           p_del_h1  = log(`0010`) - log(`1010`), 
           p_del_h2  = log(`1000`) - log(`1010`),
           p_inv_hom = log(`0101`) - log(`1010`), 
           p_inv_h1  = log(`0110`) - log(`1010`), 
           p_inv_h2  = log(`1001`) - log(`1010`),
           p_dup_hom = log(`2020`) - log(`1010`), 
           p_dup_h1  = log(`2010`) - log(`1010`), 
           p_dup_h2  = log(`1020`) - log(`1010`)) ]


# keep only entries with an SV call
# log likelihood ratio >= 1
LLR = 1
e = melt(d, 
         id.vars = c("chrom","start","end","sample","cell"), 
         measure.vars = c("p_del_h1", "p_del_h2", "p_inv_hom", "p_inv_h1", "p_inv_h2", "p_dup_hom", "p_dup_h1", "p_dup_h2"),
         variable.name = "SV_class",
         value.name    = "loglikratio",
         variable.factor = F)
e = e[loglikratio>= LLR, .SD[order(loglikratio, decreasing = T)][1,], by = .(chrom, start, end, sample, cell)]
e[, SV_class := substr(SV_class,3,nchar(SV_class))]
write.table(e, file = snakemake@output[[1]], quote=F, col.names = T, row.names = F, sep = "\t")
