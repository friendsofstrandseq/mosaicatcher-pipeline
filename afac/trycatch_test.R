library(data.table)
library(dplyr)
library(assertthat)
f_segments <- "/g/korbel2/weber/MosaiCatcher_output/POOLING/POOLING2_190822_200KB/segmentation/HGSVCxpool2/Selection_jointseg.txt"
f_counts = "/g/korbel2/weber/MosaiCatcher_output/POOLING/POOLING2_190822_200KB/counts/HGSVCxpool2/HGSVCxpool2.txt.gz"
f_segments = "/g/korbel2/weber/MosaiCatcher_files/RPE_SAMPLES/FASTQ/segmentation/RPE1-WT/Selection_jointseg.txt"
f_counts = "/g/korbel2/weber/MosaiCatcher_files/RPE_SAMPLES/FASTQ/counts/RPE1-WT/RPE1-WT.txt.gz"
chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")
CHROM <- "chr1"
### Check counts table
message(" * Reading count data ", f_counts, "...")
if (grepl("\\.gz$", f_counts)) {
    counts <- fread(f_counts)
}

# FIXME : tmp
# print(counts)
counts <- counts[counts$chrom %in% chroms, ]
# print(counts)

assert_that(
    "chrom" %in% colnames(counts),
    "start" %in% colnames(counts),
    "end" %in% colnames(counts),
    "class" %in% colnames(counts),
    "sample" %in% colnames(counts),
    "cell" %in% colnames(counts),
    "w" %in% colnames(counts),
    "c" %in% colnames(counts)
) %>% invisible()
counts[, sample_cell := paste(sample, "-", cell)]
setkey(counts, chrom, sample_cell)
bins <- unique(counts[, .(chrom, start, end)])

### Check CHROM:
assert_that(CHROM %in% unique(counts$chrom)) %>% invisible()
counts <- counts[chrom == CHROM]

print(f_segments)
message(" * Reading segmentation file from ", f_segments, "...")


# FIXME : tmp
seg <- tryCatch(
    {
        if (!is.null(f_segments)) {
            message(" * Reading segmentation file from ", f_segments, "...")
            seg <- fread(f_segments)
            seg_max <- seg[, max(bps), by = chrom]

            assert_that(
                "chrom" %in% colnames(seg),
                "bps" %in% colnames(seg)
            ) %>% invisible()
            if ("k" %in% colnames(seg)) {
                seg[, assert_that(length(unique(k)) == 1), by = .(chrom)] %>% invisible()
            }

            seg <-
                merge(seg, bins[, .N, by = chrom][, .(chrom, N = c(0, cumsum(N))[1:(.N - 1)])], by = "chrom")
            print(bins)
            print(seg)

            seg[, `:=`(from = c(1, bps[1:(.N - 1)] + 1), to = bps), by = chrom]

            seg[, `:=`(
                start = bins[from + N]$start,
                end = bins[to + N]$end
            )]

            seg[, SV_class := rep(c("bg1", "bg2"), .N)[1:.N], by = chrom]

            seg <- seg[chrom == CHROM]
        }
    },
    error = function(cond) {
        message("\n========Segmentation file processing error========\n")
        message(cond)
        message("\n\nContinuing plotting ...\n")
        return(NA)
    }
)

print(seg)
