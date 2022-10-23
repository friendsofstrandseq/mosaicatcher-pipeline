# fetch arguments
args = commandArgs(trailingOnly = T)

# # checking arguments
# if (length(args) < 3) {
#   message("Usage: Rscript GC_correction.R count-file.txt.gz gc-matrix.txt output.txt.gz optional(plots.pdf)")
#   stop()
# }
# for (a in seq(2)) {
#   if (!file.exists(args[a])) {
#     message(paste(args[a], "does not exists or cannot be found."))
#     stop()
#   }
# }

#args[1] <- '/data/r-workspace/strandseq_utils/H3JNHAFX3_MNIphotoconverted_200000_fixed.norm.txt.gz'
#args[2] <- '/data/r-workspace/strandseq_utils/GC_matrix_200000.txt'
print(snakemake@params[["gc_matrix"]])
# open files
counts <- data.table::fread(snakemake@input[["counts_vst"]], header = T)
# counts <- data.table::fread(args[1], header = T)
GC_matrix <- data.table::fread(snakemake@params[["gc_matrix"]], header = T)
# GC_matrix <- data.table::fread(args[2], header = T)
save_path <- snakemake@output[["counts_vst_gc"]]
# save_path <- args[3]

# check GC plots
# plot <- ifelse(is.na(args[4]), FALSE, TRUE)
plot <- FALSE
# if (plot) {
#   # import libraries
#   library(ggplot2)
# }

# check files
if (!all(c('cell', 'chrom', 'start', 'end', 'w', 'c') %in% colnames(counts))) {
  message("count file does not contain required columns: 'cell', 'chrom', 'start', 'end', 'w', 'c'")
  message("Usage: Rscript GC_correction.R count-file.txt.gz gc-matrix.txt output.txt.gz")
  stop()
}
if (!all(c('chrom', 'start', 'end', 'GC%') %in% colnames(GC_matrix))) {
  message("GC_matrix file does not contain required columns: 'chrom', 'start', 'end', 'GC%'")
  message("Usage: Rscript GC_correction.R count-file.txt.gz gc-matrix.txt output.txt.gz")
  stop()
}
if (! (all(unique(counts$chrom) %in% unique(GC_matrix$chrom)) &
       all(unique(counts$start) %in% unique(GC_matrix$start)) &
       all(unique(counts$end) %in% unique(GC_matrix$end)))) {
  message("bin features ('crhom', 'start', 'end') do not match between count file and GC matrix")
  message("make sure to choose files with identical bin sizes")
}


# green light message
message(paste("\ncount file:", args[1]))
message(paste("GC matrix file:", args[2]))
message(paste("savepath:", args[3]))
message("preprocessing...\n")


#################
# Preprocessing #
#################

# force cell column to factor
counts$cell <- as.factor(counts$cell)

# convert strandseq count file to count matrix
counts$tot_count <- counts$c + counts$w
count_matrix <- reshape2::dcast(counts, chrom+start+end ~ cell, value.var = "tot_count")

# add GC fractions to count_matrix not to lose matching order
count_matrix <- merge(count_matrix, GC_matrix[,c('chrom', 'start', 'end', 'GC%')], by=c('chrom', 'start', 'end'), all.x = T)
count_matrix$`GC%` <- as.numeric(count_matrix$`GC%`)

# filter unavailable GC fractions
fil <- apply(count_matrix, MARGIN = 1, FUN = (function(x) any(is.na(x))))
count_matrix.filt <- count_matrix[!fil,]
message(paste((dim(count_matrix)[1]-dim(count_matrix.filt)[1]), 'bins out of', dim(count_matrix.filt)[1], 'with NA values'))


######################
# GC bias correction #
######################

# set up the matrix
count_matrix.corrected <- data.frame(count_matrix.filt)

chrom_ids <- gsub('chr', '', unique(count_matrix.corrected$chrom))
autosomes <- !sapply(as.numeric(chrom_ids), is.na)
chr_order <- c(order(as.numeric(chrom_ids[autosomes])), sum(autosomes) + order(chrom_ids[!autosomes]))

count_matrix.corrected$chrom <- factor(count_matrix.corrected$chrom, levels=paste('chr', chrom_ids[chr_order], sep=''))

# collect plots
plots_meta <- data.frame()
plots <- list()
n <- 1

# correction chromosome by chromosome
for (chr in unique(count_matrix.corrected$chrom)) {

  # subset one chromosome
  message(paste('correction', chr))
  count_matrix.sub <- count_matrix.corrected[count_matrix.corrected$chrom == chr,]
  
  # correction cell by cell
  for (cell_n in seq(4,(ncol(count_matrix.sub)-1))) {
    
    # get counts and GC percentage for the current cell
    cell_id <- colnames(count_matrix.sub)[cell_n]
    cell_counts_raw <- as.numeric(count_matrix.sub[,(cell_id)])
    GC.frac_raw <- as.numeric(count_matrix.sub[,"GC."])
    
    # filter out count 0 for better fitting
    df <- data.frame(cbind(cell_counts_raw, GC.frac_raw))
    df <- df[df$cell_counts_raw > 0 & df$GC.frac_raw >.3, ]
    cell_counts <- df$cell_counts_raw
    GC.frac <- df$GC.frac_raw
    
    # logarithm of the counts centered on the mean
    log.norm.counts_raw <- log2( (cell_counts_raw+1) / mean(cell_counts_raw+1) )
    log.norm.counts <- log2( (cell_counts+1) / mean(cell_counts+1) )
    
    # check bins
    if (length(GC.frac) < 10 | length(log.norm.counts) < 10) {
      message(paste('skipping', chr, 'number of bins < 10'))
      next
    } 
    
    # fit curve before correction
    fit_before <- smooth.spline(GC.frac, log.norm.counts, df = 6)
    if(plot) {
      #plot_title <- 'before correction'
      #plot(GC.frac, log.norm.counts, main=plot_title, xlab="GC%", ylab="log2(counts/mean)")
      #lines(fit_before, col="red", lwd=2)
      plotdf <- data.frame(log.norm.counts, GC.frac)
      plotline <- data.frame(fit_before$x, fit_before$y)
      
      p <- ggplot(plotdf, aes(GC.frac, log.norm.counts)) +
        geom_point(size=1, alpha=.2) +
        ggtitle(paste(chr, "before")) +
        geom_line(data = plotline, aes(fit_before.x, fit_before.y), color="red")
      plots_meta[n, 1:4] <- list("index" = n, "chrom" = chr, "cell" = cell_id, "type" = "before")
      plots[[n]] <- p
      n <- n+1
    }
    
    # correct values for raw (pre-filtering) data
    pred <- predict(fit_before, GC.frac_raw)
    log.norm.corrected <- log.norm.counts_raw - pred$y
    
    # convert log of normalized counts to counts
    counts.corrected <- (2 ^ log.norm.corrected) * mean(count_matrix.sub[,cell_id]+1) -1
    
    # fit curve after correction
    if(plot) {
      # curve fit after correction is needed only for plotting
      fit_after <- smooth.spline(GC.frac_raw, log.norm.corrected, df = 6)
      #plot_title <- 'after correction'
      #plot(GC.frac_raw, log.norm.corrected, main=plot_title, xlab="GC%", ylab="log2(counts/mean)")
      #lines(fit_after, col="red", lwd=2)
      
      plotdf <- data.frame(log.norm.corrected, GC.frac_raw)
      plotline <- data.frame(fit_after$x, fit_after$y)
      p <- ggplot(plotdf, aes(GC.frac_raw, log.norm.corrected)) +
        geom_point(size=1, alpha=.2) +
        ggtitle(paste(chr, "after")) +
        geom_line(data = plotline, aes(fit_after.x, fit_after.y), color="red")
      plots_meta[n, 1:4] <- list("index" = n, "chrom" = chr, "cell" = cell_id, "type" = "after")
      plots[[n]] <- p
      n <- n+1
    }
    
    # write data to the table
    count_matrix.corrected[count_matrix.corrected$chrom == chr, cell_id] <- counts.corrected
  }
}


###################
# POST PROCESSING #
###################

# table of corrected data
counts.corrected <- reshape2::melt(count_matrix.corrected, id.vars = c("chrom", "start", "end"),
                                              measure.vars = colnames(count_matrix.corrected)[seq(4,(ncol(count_matrix.corrected)-1))], 
                                              value.name = "GC_tot_count", variable.name='cell')

# force cell column to factor
counts$cell <- as.factor(counts$cell)
counts.corrected$cell <- as.factor(counts.corrected$cell)

# merge corrected data to original data
counts.full <- merge(counts, counts.corrected, by=c("chrom", "start", "end", "cell"), all.x=T)

# ratio of the corrected and raw counts
counts.full$ratio <- counts.full$GC_tot_count/counts.full$tot_count
counts.full[is.na(counts.full$ratio)]$ratio <- 1
counts.full[is.infinite(counts.full$ratio)]$ratio <- 0

# adjust watson, crick and tot_counts
counts.full$w <- counts.full$w * counts.full$ratio
counts.full$c <- counts.full$c * counts.full$ratio
counts.full$tot_count <- counts.full$w + counts.full$c

# saving
message("saving...\n")
data.table::fwrite(counts.full[,c('chrom', 'start', 'end', 'sample', 'cell', 'w', 'c', 'class', 'tot_count')], save_path)


################
# SAVING PLOTS #
################

if (plot) {
  
  # define layout matrix
  layout_matrix <- matrix(seq(1,48), 6,8)
  # last extra squares for the dummy
  layout_matrix[5:6,8] <- 47
  
  # collect plots in right order
  n <-1
  index_order <- c()
  for (cell_id in unique(plots_meta$cell)) {
    # add dummy plot with cell name
    plots[[length(plots) + n]] <- ggplot() + 
      annotate("text", size=4, x=0, y=0, label = cell_id, angle=-90) + 
      theme_void()
    
    # retrieve indexes
    subdf <- plots_meta[plots_meta$cell == cell_id,]
    index_order <- c(index_order, as.numeric(rownames(subdf)), length(plots))
    n <- n+1
  }
  
  message("saving plots, please wait...")
  # generate the grid
  grid <- gridExtra::marrangeGrob(grobs = plots[index_order], 
                                  ncol = 8, nrow = 6, 
                                  top="GC correction plots", 
                                  layout_matrix = layout_matrix)
  
  # save plots
  ggsave(args[4], grid, width=14, height=10)
}
