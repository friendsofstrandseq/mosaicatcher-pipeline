# Median of ratios normalization

# fetch arguments
# args <- commandArgs(trailingOnly = T)

# checking arguments
# if (length(args) != 2) {
#     message("Usage: Rscript GC_correction.R count-file.txt.gz gc-matrix.txt output.txt.gz")
#     stop()
# }
# if (!file.exists(args[1])) {
#     message(paste(args[1], "does not exists or cannot be found."))
#     stop()
# }

# args[1] <- '/data/projects/strandseq_segmentation/dataset/counts/H3JNHAFX3_MNIphotoconverted_200000_fixed.txt.gz'

# open files
# counts <- data.table::fread(args[1])
counts <- data.table::fread(snakemake@input[["counts"]], header = T)
save_path <- snakemake@output[["counts_scaled"]]
# save_path <- args[3]
# min_reads <- snakemake@params[["gc_min_reads"]]

info_raw <- data.table::fread(snakemake@input[["info_raw"]], skip = 13, header = T, sep = "\t")
# info_raw <- data.table::fread(args[2])

min_reads <- min(info_raw[info_raw$pass1 == 1, ]$good) - 1
 
#################
# Preprocessing #
#################

# force cell column to factor
counts$cell <- as.factor(counts$cell)

# add tot counts
counts$tot_count <- counts$c + counts$w


# filter cells with too low counts
counts_bycell <- as.data.frame(aggregate(counts$tot_count, by = list(Category = counts$cell), FUN = sum))
sel_cells <- counts_bycell[counts_bycell$x >= as.integer(min_reads), "Category"]
if (length(sel_cells) == 0) {
    stop(paste("there are no cells with more than", min_reads, "total reads"))
}
counts <- counts[counts$cell %in% sel_cells, ]

# convert to bin matrix
count_matrix_raw <- reshape2::dcast(counts, chrom + start + end ~ cell, value.var = "tot_count")
count_matrix <- as.data.frame(count_matrix_raw)

#################
# Normalization #
#################

message(paste("library size normalization for", snakemake@input[["counts"]]))
# take the log of counts
count_matrix[4:ncol(count_matrix)] <- log(count_matrix[4:ncol(count_matrix)])

# calculate mean of the log counts per bin
count_matrix$mean_log_count <- apply(count_matrix[4:ncol(count_matrix)], MARGIN = 1, FUN = mean, na.rm = FALSE)

# filter out infinity
count_matrix <- count_matrix[!is.infinite(count_matrix$mean_log_count), ]
if (dim(count_matrix)[[1]] == 0) {
    stop("there are no common non-zero bins available across all cells.")
}

# log of counts over mean per bin
count_matrix[4:ncol(count_matrix)] <- count_matrix[4:ncol(count_matrix)] - count_matrix$mean_log_count

# median per cell of the per log of counts/mean per bin is the scaling factor
scaling_factors <- apply(count_matrix[4:ncol(count_matrix)], MARGIN = 2, FUN = median)

# raise e to scaling factor
scaling_factors <- exp(scaling_factors)

# rescale libraries
scaled_matrix <- data.table::data.table(count_matrix_raw)

for (i in colnames(scaled_matrix)[4:ncol(scaled_matrix)]) {
    scaled_matrix[[i]] <- scaled_matrix[[i]] / scaling_factors[i]
}

##################
# Postprocessing #
##################

# wide to long
norm_tot_counts <- reshape2::melt(scaled_matrix,
    id.vars = c("chrom", "start", "end"),
    measure.vars = colnames(scaled_matrix)[4:ncol(scaled_matrix)],
    value.name = "norm_tot_count", variable.name = "cell"
)

# merge with counts
cols <- colnames(counts)
counts <- merge(counts, norm_tot_counts, by = c("chrom", "start", "end", "cell"), all.x = TRUE)

# adjust W and C counts to normalized counts
counts$ratio <- counts$norm_tot_count / counts$tot_count
counts$ratio[is.na(counts$ratio)] <- 0
counts$w <- counts$w * counts$ratio
counts$c <- counts$c * counts$ratio
counts$tot_count <- counts$norm_tot_count
# fill na
counts$w[is.na(counts$w)] <- 0
counts$c[is.na(counts$c)] <- 0
counts$tot_count[is.na(counts$tot_count)] <- 0


# saving
message("saving...\n")
data.table::fwrite(counts[, ..cols], save_path)
