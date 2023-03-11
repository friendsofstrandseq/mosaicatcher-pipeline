args <- commandArgs(trailingOnly = TRUE)


if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update = FALSE)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", update = FALSE)
quit(save = "no")
