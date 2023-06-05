# log <- file(snakemake@log[[1]], open = "wt")
# log <- file("test.log", open = "wt")
# sink(file = log, type = .snakemake/conda/822b077b26285271f91dfda88e8c86b5"output")

# IMPORTS

library(ComplexHeatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = T)

data_file <- args[1]
output_file <- args[2]
stringent_lenient <- args[3]
filter <- args[4]
sample_name <- args[5]

print(data_file)
print(output_file)
print(stringent_lenient)
print(filter)

# pdf("TEST_R_dev.pdf", width = 20, height = 10)
pdf(output_file, width = 20, height = 10)

# Chromosome order
chrOrder <-
    c(paste("chr", 1:22, sep = ""), "chrX")

# Load SV data

# data_file = "../stringent_filterTRUE.tsv"
# data_file <- "/scratch/tweber/SCO_COURSE/HJ_MIXTURE_RPE1_Mix/RPE1_Mix/mosaiclassifier/sv_calls/stringent_filterTRUE.tsv"
# data1 <- read.table("../lenient_filterFALSE.tsv",
data1 <- read.table(data_file,
    sep = "\t",
    header = T,
    comment.char = ""
)
# head(data1)

# Create Dataframe for chromosomes missing SVs

chrom <- as.vector(setdiff(chrOrder, data1$chrom))
start <- rep(0, length(setdiff(chrOrder, data1$chrom)))
end <- rep(0, length(setdiff(chrOrder, data1$chrom)))
sample <- rep(data1$sample[1][1], length(setdiff(chrOrder, data1$chrom)))
cell <- rep(data1$cell[1][1], length(setdiff(chrOrder, data1$chrom)))
class <- rep("NA", length(setdiff(chrOrder, data1$chrom)))
scalar <- rep(0, length(setdiff(chrOrder, data1$chrom)))
num_bins <- rep(0, length(setdiff(chrOrder, data1$chrom)))
sv_call_name <- rep("none", length(setdiff(chrOrder, data1$chrom)))
sv_call_haplotype <- rep(0, length(setdiff(chrOrder, data1$chrom)))
sv_call_name_2nd <- rep("NA", length(setdiff(chrOrder, data1$chrom)))
sv_call_haplotype_2nd <- rep(0, length(setdiff(chrOrder, data1$chrom)))
llr_to_ref <- rep(0, length(setdiff(chrOrder, data1$chrom)))
llr_to_2nd <- rep(0, length(setdiff(chrOrder, data1$chrom)))
af <- rep(0, length(setdiff(chrOrder, data1$chrom)))

data1_missing <- data.frame(
    chrom,
    start,
    end,
    sample,
    cell,
    class,
    scalar,
    num_bins,
    sv_call_name,
    sv_call_haplotype,
    sv_call_name_2nd,
    sv_call_haplotype_2nd,
    llr_to_ref,
    llr_to_2nd,
    af
)


# Bind existing dataframe and new one

data1 <- rbind(data1, data1_missing)

data1$chrom <-
    factor(data1$chr, levels = chrOrder)


data1 <- data1[order(data1$chrom), ]
data1$pos <- paste0(data1$chrom, "_", data1$start, "_", data1$end)

# Select subset of the dataframe
lite_data <- select(data1, c("pos", "cell", "sv_call_name"))

# Get colors / chrom
set.seed(2)
n <- length(unique(data1$chrom))
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1, n), col = sample(col_vector, n))

chrom <- unique(data1$chrom)
colors_chroms <- sample(col_vector, length(chrom))
# Instanciate data$color column
data1$color <- "NA"
# Iterate over chrom to attribute color
for (i in 1:length(chrom)) {
    chrom_index_list <- which(data1$chrom == chrom[i])
    data1[chrom_index_list, "color"] <- colors_chroms[i]
}
dd <- unique(select(data1, c("chrom", "color")))

# range01 <- function(x) {
#     100 * ((x - min(x)) / (max(x) - min(x)))
# }

## LLR & CLUSTERING

# Create subset for clustering

lite_data_clustering <- select(data1, c("pos", "cell", "llr_to_ref"))
lite_data_clustering[c("llr_to_ref")][sapply(lite_data_clustering[c("llr_to_ref")], is.infinite)] <- max(lite_data_clustering$llr_to_ref[is.finite(lite_data_clustering$llr_to_ref)])
lite_data_clustering[is.na(lite_data_clustering)] <- 0
# lite_data_clustering$llr_to_ref <- range01(lite_data_clustering$llr_to_ref)

# Pivot dataframe into matrix
lite_data_pivot_clustering <- lite_data_clustering %>%
    pivot_wider(
        names_from = "cell",
        values_from = "llr_to_ref"
    )
# Transpose
t_lite_data_pivot_clustering <- t(lite_data_pivot_clustering)

colnames(t_lite_data_pivot_clustering) <- t_lite_data_pivot_clustering[1, ]
t_lite_data_pivot_clustering <- t_lite_data_pivot_clustering[-1, ]
t_lite_data_pivot_clustering[is.na(t_lite_data_pivot_clustering)] <- 0

# Turn into numeric matrix
t_lite_data_pivot_clustering_num <- matrix(as.double(t_lite_data_pivot_clustering), ncol = ncol(t_lite_data_pivot_clustering))
rownames(t_lite_data_pivot_clustering_num) <- rownames(t_lite_data_pivot_clustering)
colnames(t_lite_data_pivot_clustering_num) <- colnames(t_lite_data_pivot_clustering)

# Plot options
options(repr.plot.width = 20, repr.plot.height = 12)

anno_colors <- list(Chroms = unique(data1$chrom))

col_annotation <- sapply(strsplit(lite_data_pivot_clustering$pos, "_"), `[`, 1)

col_test <- factor(sapply(strsplit(colnames(t_lite_data_pivot_clustering_num), "_"), `[`, 1), levels = unique(sapply(strsplit(colnames(t_lite_data_pivot_clustering_num), "_"), `[`, 1)))

# print(paste0("Chromosome size unscaled LLR heatmap (Sample : ", sample, ", Methods used: ", stringent_lenient, ", Filter used: ", filter, ")"))

# print(t_lite_data_pivot_clustering_num)

cl_h <- Heatmap(as.matrix(t_lite_data_pivot_clustering_num),
    name = "LLR", col = RColorBrewer::brewer.pal(name = "Reds", n = 9),
    # column_title = "a discrete numeric matrix",
    rect_gp = gpar(col = "white", lwd = 1.5),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
        foo = anno_block(gp = gpar(fill = 2:24))
    ),
    column_split = col_test,
    width = unit(32, "cm"), height = unit(20, "cm"),
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 4),
    column_title_gp = gpar(fontsize = 10),
    cluster_columns = FALSE,
    column_gap = unit(2, "mm"),
    cluster_column_slices = FALSE,
    column_title_rot = 90,
)
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")
draw(
    cl_h,
    # row_title = "Three heatmaps, row title", row_title_gp = gpar(col = "red"),
    # column_title = paste0("Chromosome size unscaled LLR heatmap (Sample : ", snakemake@wildcards[["sample"]], ", Methods used: ", snakemake@wildcards[["method"]], ", Filter used: ", snakemake@wildcards[["filter"]], ")"), column_title_gp = gpar(fontsize = 16)
    column_title = paste0("Chromosome size unscaled LLR heatmap (Sample : ", sample_name, ", Methods used: ", stringent_lenient, ", Filter used: ", filter, ")"), column_title_gp = gpar(fontsize = 16)
)

## CATEGORICAL

# Turn data into a matrix
lite_data_pivot <- lite_data %>%
    pivot_wider(
        names_from = "cell",
        values_from = "sv_call_name"
    )

# Transpose
t_lite_data_pivot <- t(lite_data_pivot)
colnames(t_lite_data_pivot) <- t_lite_data_pivot[1, ]
t_lite_data_pivot <- t_lite_data_pivot[-1, ]

# SV list
sv_list <-
    c(
        "none",
        "del_h1",
        "del_h2",
        "del_hom",
        "dup_h1",
        "dup_h2",
        "dup_hom",
        "inv_h1",
        "inv_h2",
        "inv_hom",
        "idup_h1",
        "idup_h2",
        "complex"
    )

# SV type colors
colors <-
    structure(
        c(
            "grey",
            "#77AADD",
            "#4477AA",
            "#114477",
            "#CC99BB",
            "#AA4488",
            "#771155",
            "#DDDD77",
            "#AAAA44",
            "#777711",
            "#DDAA77",
            "#AA7744",
            "#774411"
        ),
        names = sv_list
    )

# Fill NA with none
t_lite_data_pivot[is.na(t_lite_data_pivot)] <- "none"

anno_colors <- list(Chroms = unique(data1$chrom))
col_annotation <- as.data.frame(sapply(strsplit(lite_data_pivot$pos, "_"), `[`, 1))
colnames(col_annotation) <- "Chroms"
col_test <- factor(sapply(strsplit(colnames(t_lite_data_pivot), "_"), `[`, 1), levels = unique(sapply(strsplit(colnames(t_lite_data_pivot), "_"), `[`, 1)))


cat_h <- Heatmap(as.matrix(t_lite_data_pivot),
    name = "SV type", col = colors,
    # column_title = "a discrete numeric matrix",
    rect_gp = gpar(col = "white", lwd = 1.5),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
        foo = anno_block(gp = gpar(fill = 2:24))
    ),
    column_split = col_test,
    width = unit(32, "cm"), height = unit(20, "cm"),
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 4),
    column_title_gp = gpar(fontsize = 10),
    column_gap = unit(2, "mm"),
    # column_order = order(as.numeric(sapply(strsplit(gsub("chr", "", colnames(t_lite_data_pivot)), "_"), `[`, 1))),
    row_order = row_order(cl_h),
    column_title_rot = 90,
    # use_raster = TRUE, raster_by_magick = TRUE, raster_quality=10
)
ht_opt$TITLE_PADDING <- unit(c(8.5, 8.5), "points")
draw(
    cat_h,
    # row_title = "Three heatmaps, row title", row_title_gp = gpar(col = "red"),
    # column_title = paste0("Chromosome size unscaled categorical heatmap (Sample : ", snakemake@wildcards[["sample"]], ", Methods used: ", snakemake@wildcards[["method"]], ", Filter used: ", snakemake@wildcards[["filter"]], ")"), column_title_gp = gpar(fontsize = 16)
    column_title = paste0("Chromosome size unscaled categorical heatmap (Sample : ", sample_name, ", Methods used: ", stringent_lenient, ", Filter used: ", filter, ")"), column_title_gp = gpar(fontsize = 16)
)


# Export clustered row order to output in order to use it in python script
row_order <- row_order(cl_h)
cell <- rownames(t_lite_data_pivot)[row_order]
index <- seq(1, length(cell))
cluster_order_df <- data.frame(index, row_order, cell)
# write.table(cluster_order_df, file = "test.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cluster_order_df, file = "TEST_cluster_order_df.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
