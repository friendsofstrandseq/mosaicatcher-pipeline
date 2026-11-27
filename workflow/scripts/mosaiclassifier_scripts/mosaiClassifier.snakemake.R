# Send both output and messages (warnings, assert_that errors) to the rule log;
# otherwise errors only reach the SLURM log, which is purged.
log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = "output")
sink(file = log, type = "message")
library(data.table)
library(assertthat)
source("workflow/scripts/mosaiclassifier_scripts/mosaiClassifier/mosaiClassifier.R")
system("LC_CTYPE=C")

# Currently read files from the Snakemake pipeline
counts <- fread(paste("zcat", snakemake@input[["counts"]]))
info <- fread(snakemake@input[["info"]])
strand <- fread(snakemake@input[["states"]])
segs <- fread(snakemake@input[["bp"]])


# FIXME : tmp solution to fix error : Segments must covered all bins, which happen for small scaffolds
chroms <- snakemake@config[["chromosomes"]]

counts <- counts[counts$chrom %in% chroms, ]
strand <- strand[strand$chrom %in% chroms, ]
segs <- segs[segs$chrom %in% chroms, ]

# A chromosome without any segment cannot be classified: drop it from counts and
# strand states. This happens when every bin of the chromosome is class "None",
# e.g. chrY in female mouse cells (no read at all), for which `mosaicatcher
# segment` emits no breakpoint. Keeping it makes addCountsPerSegment() fail with
# "Segments should cover all bins".
seg_chroms <- unique(segs$chrom)
dropped <- setdiff(unique(counts$chrom), seg_chroms)
if (length(dropped) > 0) {
  message(
    "[MosaiClassifier] No segment for ", paste(dropped, collapse = ", "),
    ": chromosome(s) dropped"
  )
  counts <- counts[counts$chrom %in% seg_chroms, ]
  strand <- strand[strand$chrom %in% seg_chroms, ]
}

# haplotypeMode?
if ("CW" %in% strand$class) {
  haplotypeMode <- T
} else {
  haplotypeMode <- F
}

d <- mosaiClassifierPrepare(counts, info, strand, segs)
e <- mosaiClassifierCalcProbs(d, maximumCN = 4, haplotypeMode = haplotypeMode)

saveRDS(e, file = snakemake@output[[1]])