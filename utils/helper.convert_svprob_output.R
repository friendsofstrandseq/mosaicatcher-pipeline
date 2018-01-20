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

d <- d[, .(chrom = chr,
           start = format(start,scientific=F),
           end   = format(end,scientific=F),
           sample, cell, 
           type = toupper(types), w = Wcount, c = Ccount,
           p_cn0 = CN0, p_cn1 = CN1, p_cn2 = CN2, p_cn3 = CN3, p_cn4 = CN4,
           p_ref = `1010`, p_del_hom = `0000`, p_del_h1 = `0010`, p_del_h2 = `1000`,
           p_inv_hom = `0101`, p_inv_h1 = `0110`, p_inv_h2 = `1001`,
           p_dup_hom = `2020`, p_dup_h1 = `2010`, p_dup_h2 = `1020`,
           p_idup_h1 = `1110`, p_idup_h2 = `1011`)]
d
write.table(d, file = snakemake@output[[1]], quote=F, col.names = T, row.names = F, sep = "\t")
