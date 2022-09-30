library(edgeR)
library(ggplot2)

# fetch arguments
args = commandArgs(trailingOnly = T)

filter <- ifelse('filter' %in% args, TRUE, FALSE)
chosen_transform <- ifelse('anscombe' %in% args, 'anscombe', 
                           ifelse('laubschner' %in% args, 'laubschner', 'anscombe'))

# open count file
#args[1] <- '/data/strandseq_datasets/lib_test/HCLLVAFX3_iPSCs.200000.LIB.txt.gz'
# counts_raw <- data.table::fread(args[1])
counts_raw <- data.table::fread(snakemake@input[["counts"]])

# force cell column to factor
counts_raw$cell <- as.factor(counts_raw$cell)

# fuse bin coordinates
counts_raw$bin <- paste(counts_raw$chrom, counts_raw$start, counts_raw$end, sep='_')

# add tot counts
counts_raw$tot_count <- counts_raw$c + counts_raw$w

# sum of counts per cell
sum_counts <- aggregate(counts_raw$tot_count, by= list(counts_raw$cell), FUN = sum)
colnames(sum_counts) <- c('cell', 'sum')


# centromere read spikes exclusion

exclude_centromeres <- function(counts_raw, exclusion_list) {
  exclusion_list <- GenomicRanges::GRanges(seqnames = exclusion_list$chrom, ranges = IRanges::IRanges(start = exclusion_list$start, end=exclusion_list$end))
  counts_iranges <- GenomicRanges::GRanges(seqnames = counts_raw$chrom, ranges = IRanges::IRanges(start = counts_raw$start, end=counts_raw$end))

  overlaps <- data.table::as.data.table(GenomicRanges::findOverlaps(exclusion_list, counts_iranges))
  counts <- counts_raw
  counts <- counts[-overlaps$subjectHits,]
  return(counts)
}

#args[3] <- '/data/r-workspace/strandseq_utils/large_centromere_blacklist.csv'

if (!is.na(args[3])) {
  message('blacklisting...')
  counts <- exclude_centromeres(counts_raw, data.table::fread(args[3]))
} else {
  counts <- counts_raw
}

# convert to bin count matrix
to_matrix <- function(counts) {
  mat_tot <- reshape2::dcast(counts, bin ~ cell, value.var = "tot_count")
  rownames(mat_tot) <- mat_tot$bin
  mat_tot <- mat_tot[,2:ncol(mat_tot)]
  
  return(mat_tot)
}

mat_tot <- to_matrix(counts)

generate_dgelist <- function(mat_tot, filter) {
  
  # create the DGEList object
  # all samples belong to group 1
  y.raw <- DGEList(as.matrix(mat_tot), group = rep(1, ncol(mat_tot)))
  
  # filter out bins with too low counts to be informative
  if (filter) {
    keep.exprs <- edgeR::filterByExpr(y.raw, group=y.raw$samples$group)
    y <- y.raw[keep.exprs, keep.lib.sizes=FALSE]
    message(paste('filtering out', (dim(y.raw)[1]-dim(y)[1]), 'bins out of', dim(y.raw)[1], "due to low count"))
  } else {
    y <- y.raw
    message('no filtering')
  }
  
  # calculate library size and composition normalization
  y <- calcNormFactors(y)
  
  return(y)
}

estimate_dispersion <- function(y) {
  
  # estimate common dispersion
  # default settings for DGEList objects
  message('Estimating common dispersion with edgeR...')
  disp <- estimateDisp(y, design=NULL, prior.df=NULL, trend.method="locfit", tagwise=TRUE,
                       span=NULL, min.row.sum=5, grid.length=21, grid.range=c(-10,10), robust=FALSE, 
                       winsor.tail.p=c(0.05,0.1), tol=1e-06)
  phi <- disp$common.dispersion
  message(paste("phi =",phi))
  write(paste(args[1], phi, sep="\t"), './log.txt', append=TRUE)
  
  return(phi)
}

y <- generate_dgelist(mat_tot, filter)
phi <- estimate_dispersion(y)


# VARIANCE STABILIZING TRANSFORMATION

# TRANSFORMS
# for negative binomial distributions
# Anscombe, 1948
# Laubschner, 1961
anscombe_transform <- function(x, phi) {
  a <- x + (3/8)
  b <- (1/phi) - (3/4)
  c <- sqrt(a/b)
  y <- sinh(c)
  return(y)
}
laubscher_transform <- function(x, phi) {
  a <- sqrt(phi)
  b <- sinh( sqrt( x/phi ) )
  c <- sqrt( phi-1 )
  d <- anscombe_transform(x, phi)
  y <- a*b + c*d
  return(y)
}
transform_list <- list('anscombe'= anscombe_transform, 'laubscher'= laubscher_transform)
transform <- transform_list[[chosen_transform]]

message(paste("Transforming data with", chosen_transform, "VST"))
#ans <- transform(y$counts, phi)
ans <- transform(to_matrix(counts_raw), phi)

# convert to count table
ans <- as.data.frame(ans)
ans$bin <- rownames(ans)
d <- reshape2::melt(ans, id.vars='bin', measure.vars=colnames(ans), variable.name='cell', value.name = 'tot_count_corr')

sum2 <- apply(ans[,1:ncol(ans)-1], MARGIN = 2, FUN = sum)
sum2 <- data.table::data.table(cell = names(sum2), sum2 = sum2)
sum_counts <- merge(sum_counts, sum2, by = c('cell'), all.x = T)
sum_counts$f <- sum_counts$sum / sum_counts$sum2

d <- merge(d, sum_counts[,c('cell', 'f')], by = c('cell'), all.x = T)
d$tot_count_corr <- as.numeric(d$tot_count_corr) * d$f
d <- d[,c('cell', 'bin', 'tot_count_corr')]

bins <- stringr::str_split_fixed(d$bin, '_', 3)
colnames(bins) <- c('chrom', 'start', 'end')
d[c('chrom', 'start', 'end')] <- bins
d$start <- as.numeric(d$start)
d$end <- as.numeric(d$end)

merged.raw <- merge(counts_raw, d[c('chrom', 'start', 'end', 'cell', 'tot_count_corr')], 
                    by=c('chrom', 'start', 'end', 'cell'), all.x=T)

fil <- apply(merged.raw, MARGIN = 1, FUN = (function(x) any(is.na(x))))
#merged.raw[fil, 'tot_count_corr'] <- 0
merged <- data.table::data.table(merged.raw)

# adjust W and C counts to corrected tot_counts
merged$tot_count_corr <- as.numeric(merged$tot_count_corr)
merged$ratio <- merged$tot_count_corr / merged$tot_count
merged$w <- merged$w * merged$ratio
merged$c <- merged$c * merged$ratio
merged$w[is.na(merged$w)] <- 0
merged$c[is.na(merged$c)] <- 0
merged$tot_count <- merged$tot_count_corr

message('saving...')
df <- data.table::data.table(merged[,c('chrom', 'start', 'end', 'sample', 'cell', 'w', 'c', 'class', 'tot_count')])

#s1 <- aggregate(counts_raw$tot_count, by= list(counts_raw$cell), FUN = sum)
#s2 <- aggregate(df$tot_count, by= list(df$cell), FUN = sum)
#s3 <- merge(s1, s2, by=c('Group.1'), all.x=T)
#s3$d <- s3['x.x'] - s3['x.y']
#df

#args[2] <- '/data/r-workspace/strandseq_utils/vst_test.txt.gz'
# data.table::fwrite(df, args[2])
data.table::fwrite(df, snakemake@output[["counts_vst"]])

#scTRIPmultiplot::multiplot(args[2], chromosome = 'chr8', cell_id = 'MNIx260521x01PE20505')
#scTRIPmultiplot::multiplot(args[1], chromosome = 'chr8', cell_id = 'MNIx260521x01PE20505')
#path <- '/data/strandseq_datasets/H22WGAFX3_MNIx260521.200000.VST.txt.gz'
#scTRIPmultiplot::multiplot(path, chromosome = 'chr8', cell_id = 'MNIx260521x01PE20505')
