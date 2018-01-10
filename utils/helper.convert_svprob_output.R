sink(snakemake@log[[1]])
library(data.table)
library(assertthat)

f_maryam = snakemake@input[["probs"]]
f_maryam
f_info   = snakemake@input[["info"]]
f_info
sample   = snakemake@params[["sample_name"]]
sample

d = fread(f_maryam, header = T)
d
f = fread(f_info)
f

# Map cell number to cell name
assert_that(max(d$cells) <= nrow(f))
d$cell = f$cell[d$cells]
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
