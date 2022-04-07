# sink(snakemake@log[[1]])
library(data.table)
library(assertthat)
source("utils/mosaiClassifier/mosaiClassifier.R")


# Currently read files from the Snakemake pipeline
counts <- fread(paste("zcat", "/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG/counts/RPE1-WT/100000.txt.gz"))
info <- fread("/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG/counts/RPE1-WT/100000.info")
strand <- fread("/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG/strand_states/RPE1-WT/100000.selected_j0.1_s0.5_scedist20/final.txt")
segs <- fread("/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG/segmentation2/RPE1-WT/100000.selected_j0.1_s0.5_scedist20.txt")

chroms <- c("chr1", "chr2")

counts = counts[counts$chrom %in% chroms, ]
strand = strand[strand$chrom %in% chroms, ]
segs = segs[segs$chrom %in% chroms, ]

print(counts)
print(strand)
print(segs)

# haplotypeMode?
if ("CW" %in% strand$class) {
    haplotypeMode <- T
} else {
    haplotypeMode <- F
}

d <- mosaiClassifierPrepare(counts, info, strand, segs)
print(d)
e <- mosaiClassifierCalcProbs(d, maximumCN = 4, haplotypeMode = haplotypeMode)

saveRDS(e, file = "/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG/sv_probabilities/RPE1-WT/100000.selected_j0.1_s0.5_scedist20/probabilities.Rdata")