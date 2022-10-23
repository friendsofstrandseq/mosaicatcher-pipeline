# fetch arguments
args <- commandArgs(trailingOnly = T)

# args[1] <- '/data/r-workspace/strandseq_utils/test_data/H3JNHAFX3_MNIphotoconverted.200000.VST.GC.txt.gz'

# open counts
# counts <- data.table::fread(args[1])
counts <- data.table::fread(snakemake@input[["counts_vst_gc"]])


# calculate scale factor essentially checks the order of magnitude (log10)
# and rescales to the 100s
scale_factor <- (floor(log10(median(counts$w + counts$c))) - 2) * -1

message(paste("scale factor:", scale_factor))

# scaling
counts$w <- as.integer(counts$w * 10^scale_factor)
counts$c <- as.integer(counts$c * 10^scale_factor)

if ("tot_count" %in% colnames(counts)) {
    counts$tot_count <- as.integer(counts$w + counts$c)
}

# filename
# tokens <- strsplit(args[1], "/")[[1]]
# basename <- tokens[length(tokens)]
# basename <- strsplit(args[1], ".txt.gz")[[1]]

# save
# data.table::fwrite(counts, args[2])
data.table::fwrite(counts, snakemake@output[["counts_vst_gc_scaled"]])
