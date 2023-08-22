# set arguments
input_path <- snakemake@input[["counts_scaled_gc"]]
save_path <- snakemake@output[["counts_scaled_gc_vst"]]
chosen_transform <- "anscombe"
plot <- TRUE
rescale <- TRUE

# PRE-PROCESSING DATA
# open count file
counts_raw <- data.table::fread(input_path)

# force cell column to factor
counts_raw$cell <- as.factor(counts_raw$cell)

# fuse bin coordinates
counts_raw$bin <- paste(counts_raw$chrom, counts_raw$start, counts_raw$end, sep = "_")

# add tot counts
counts_raw$tot_count <- counts_raw$c + counts_raw$w

counts <- counts_raw

# convert to bin count matrix
to_matrix <- function(counts) {
  # fuse bin coordinates
  counts$bin <- paste(counts$chrom, counts$start, counts$end, sep = "_")
  mat_tot <- reshape2::dcast(counts, bin ~ cell, value.var = "tot_count")
  rownames(mat_tot) <- mat_tot$bin
  mat_tot <- mat_tot[, 2:ncol(mat_tot)]

  return(mat_tot)
}

# VARIANCE STABILIZING TRANSFORMATION

# TRANSFORMS
# for negative binomial distributions
# Anscombe, 1948
# Laubschner, 1961
anscombe_transform <- function(x, phi) {
  a <- x + (3 / 8)
  b <- (1 / phi) - (3 / 4)
  c <- sqrt(a / b)
  y <- asinh(c)
  return(y)
}
laubscher_transform <- function(x, phi) {
  a <- sqrt(phi)
  b <- asinh(sqrt(x / phi))
  c <- sqrt(phi - 1)
  d <- anscombe_transform(x, phi)
  y <- a * b + c * d
  return(y)
}
transform_data <- function(counts, transform, phi) {
  cols <- colnames(counts)
  counts$tot_count_corr <- transform(counts$tot_count, phi)
  counts$f <- counts$tot_count_corr / counts$tot_count
  counts$w <- counts$w * counts$f
  counts$c <- counts$c * counts$f
  counts$tot_count <- counts$tot_count_corr
  counts <- counts[, ..cols]

  return(counts)
}
transform_list <- list("anscombe" = anscombe_transform, "laubscher" = laubscher_transform)
transform <- transform_list[[chosen_transform]]


disp_score <- function(counts, transform, phi, design = NULL) {
  counts$tot_count <- transform(counts$tot_count, phi)
  mat <- to_matrix(counts)
  # if multiple samples are present design matrix can be used
  if (is.null(design)) {
    design <- matrix(1, ncol = 1, nrow = ncol(mat))
  }
  res <- as.matrix(mat) %*% MASS::Null(design)
  rsd <- sqrt(rowMeans(res * res))
  score <- sd(rsd) / mean(rsd)
  return(score)
}

message(paste("Transforming data with", chosen_transform, "VST"))

# estimate dispersion by residual variance

opt <- optimize(disp_score, counts = counts, transform = transform, interval = c(0.00001, 1))
phi <- opt$minimum
message(paste("Estimated dispersion - phi: ", phi))

# correction
corr_counts <- transform_data(counts, transform, phi)
corr_counts <- data.table::data.table(corr_counts[, c("chrom", "start", "end", "sample", "cell", "w", "c", "class", "tot_count")])

rescale_data <- function(counts_original, counts_transformed) {
  rescaled <- counts_transformed
  rescaled_med <- aggregate(rescaled$tot_count, list(rescaled$cell), FUN = median)
  original_med <- aggregate(counts_original$tot_count, list(counts_original$cell), FUN = median)
  m <- merge(x = original_med, y = rescaled_med, by = "Group.1", suffixes = c("_raw", "_norm"))
  m[["f"]] <- m[["x_raw"]] / m[["x_norm"]]

  rescaled <- merge(rescaled, m[c("Group.1", "f")], by.x = "cell", by.y = "Group.1")

  rescaled$tot_count <- rescaled$tot_count * rescaled$f
  rescaled$w <- rescaled$w * rescaled$f
  rescaled$c <- rescaled$c * rescaled$f
  return(rescaled)

}

if (rescale == TRUE) {
  corr_counts <- rescale_data(counts_raw, corr_counts)
} 


message("saving...")
data.table::fwrite(corr_counts, save_path)

if (plot) {
  library(ggplot2)
  library(ggpubr)

  merge_bins <- function(df, bin_size = 3e6) {
    df <- df[with(df, order(cell, chrom, start))]

    df$bin_group <- df$start %/% bin_size


    w <- aggregate(df$w, by = list(df$cell, df$chrom, df$bin_group), FUN = sum)
    c <- aggregate(df$c, by = list(df$cell, df$chrom, df$bin_group), FUN = sum)
    s <- aggregate(df$start, by = list(df$cell, df$chrom, df$bin_group), FUN = function(x) x[[1]])
    e <- aggregate(df$end, by = list(df$cell, df$chrom, df$bin_group), FUN = function(x) x[[length(x)]])
    cl <- aggregate(df$class, by = list(df$cell, df$chrom, df$bin_group), FUN = function(x) names(sort(table(x), decreasing = TRUE))[[1]])

    m <- merge(w, c, by = c("Group.1", "Group.2", "Group.3"))
    m <- merge(m, s, by = c("Group.1", "Group.2", "Group.3"))
    m <- merge(m, e, by = c("Group.1", "Group.2", "Group.3"))
    m <- merge(m, cl, by = c("Group.1", "Group.2", "Group.3"))
    colnames(m) <- c("cell", "chrom", "bin_group", "w", "c", "start", "end", "class")
    m$sample <- unique(df$sample)[[1]]
    m <- data.table::as.data.table(m)
    m <- m[with(m, order(cell, chrom, start))]

    m$tot_count <- m$w + m$c
    return(m)
  }

  wf_plot <- function(m) {
    m$wf <- m$w / m$tot_count

    p1 <- ggplot(m, aes(x = tot_count, y = wf)) +
      geom_point(size = 1, alpha = .1, shape = 16) +
      xlab("tot count") +
      ylab("watson fraction") +
      ylim(0, 1)
    return(p1)
  }

  p1 <- ggplot(counts_raw, aes(x = tot_count)) +
    geom_histogram(bins = 256) +
    ggtitle("raw") +
    xlab("read count") +
    ylab("bin count")

  p2 <- ggplot(corr_counts, aes(x = tot_count)) +
    geom_histogram(bins = 256) +
    ggtitle(paste(chosen_transform, "VST")) +
    xlab("read count") +
    ylab("bin count")
  
  m <- merge_bins(counts_raw)
  p3 <- wf_plot(m) + ggtitle('raw')
  
  n <- merge_bins(corr_counts)
  p4 <- wf_plot(n) + ggtitle(paste(chosen_transform, "VST"))



  m <- merge_bins(counts_raw)
  p3 <- wf_plot(m) + ggtitle("raw")

  n <- merge_bins(corr_counts)
  p4 <- wf_plot(n) + ggtitle(paste(chosen_transform, "VST"))


  corr_plot <- ggarrange(p1, p2, p3, p4)

  ggsave(snakemake@output[["plot"]], corr_plot, width = 12, height = 6)
}