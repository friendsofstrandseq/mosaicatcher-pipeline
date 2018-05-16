library(dplyr)
library(data.table)
library(assertthat)

# Prepare function to assess the strand_state of a cell/chrom/region
addStrandStatesInner <- function(sample_, cell_, chrom_, from_, to_) {
  x = strand[sample == sample_ & cell == cell_ & chrom == chrom_]
  assert_that(nrow(x)>0)
  if (nrow(x) == 1) return (x$class)
  min_pos = bins[chrom == chrom_]$start[from_]
  max_pos = bins[chrom == chrom_]$end[to_]
  x = x[start <= min_pos & end >= max_pos]
  if (nrow(x) == 1) return (x$class)
  return ("sce")
}

addStrandStates <- function(counts, states) {
  bins <- unique(counts[, .(chrom, start, end)])
  @ work here
  
}