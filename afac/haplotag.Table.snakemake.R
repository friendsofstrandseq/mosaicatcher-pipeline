
library(devtools)
load_all("/g/korbel2/weber/Gits/Rsamtools/")
load_all("/g/korbel2/weber/Gits/GenomicAlignments/")

source("utils/haplotagTable.R")


# tab <- getHaplotagTable2(
#     bedFile = "/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG_chr21/haplotag/bed/RPE1-WT/100000.selected_j0.1_s0.5_scedist20.bed",
#     bam.file = "/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG_chr21/haplotag/bam/RPE1-WT/100000.selected_j0.1_s0.5_scedist20/RPE1WTPE20490.sort.mdup.bam",
#     file.destination = "/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG_chr21/test.tsv"
# )

tab <- getHaplotagTable2(
    bedFile = "/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG/haplotag/bed/RPE1-WT/100000.selected_j0.1_s0.5_scedist20.bed",
    bam.file = "/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG/haplotag/bam/RPE1-WT/100000.selected_j0.1_s0.5_scedist20/RPE1WTPE20490.sort.mdup.bam",
    file.destination = "/g/korbel2/weber/MosaiCatcher_output/MosaiCatcher_output_sample_KG/test.tsv"
)