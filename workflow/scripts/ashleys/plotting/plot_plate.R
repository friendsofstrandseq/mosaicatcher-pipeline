library(platetools)
library(ggplot2)
library(viridis)
library(dplyr)

## collect ASHLEYS prediction and count files
# ashleys_data <- read.table(file = "/scratch/tweber/DATA/TMP/labels384.tsv", sep = "\t", header = TRUE)
ashleys_data <- read.table(file = snakemake@input[["labels"]], sep = "\t", header = TRUE)
plate_type <- nrow((ashleys_data))
ashleys_data <- dplyr::arrange(ashleys_data, cell)
colnames(ashleys_data)[1] <- "ashleys_id"

Well_position <- character()

if (plate_type == 96) {
    for (i in 1:12)
    {
        for (j in 1:8)
        {
            tmp <- paste0(LETTERS[j], i)
            Well_position <- c(Well_position, tmp)
        }
    }
} else if (plate_type == 384) {
    for (i in 1:24)
    {
        for (j in 1:16)
        {
            tmp <- paste0(LETTERS[j], i)
            Well_position <- c(Well_position, tmp)
        }
    }
}

# pdf("TEST_ashleys_plate_predictions.pdf")
pdf(snakemake@output[["predictions"]])

ashleys_data$Well_position <- Well_position

raw_map(
    data = ashleys_data$prediction,
    well = ashleys_data$Well_position,
    plate = plate_type
) +
    scale_fill_distiller(type = "div", palette = "RdYlGn", direction = 1) +
    # ggtitle(paste0("Sample: TEST | ASHLEYS binary predictions (cutoff=0.5)"))
    ggtitle(paste0("Sample: ", snakemake@wildcards[["sample"]], " | ASHLEYS binary predictions (cutoff=", snakemake@config[["ashleys_threshold"]], ")"))

dev.off()

# pdf("TEST_ashleys_plate_probabilities.pdf")
pdf(snakemake@output[["probabilities"]])
ashleys_data$Well_position <- Well_position

raw_map(
    data = ashleys_data$probability,
    well = ashleys_data$Well_position,
    plate = plate_type
) +
    scale_fill_distiller(type = "div", palette = "RdYlGn", direction = 1) +
    # ggtitle(paste0("Sample: ", "TEST", " | ASHLEYS probabilities"))
    ggtitle(paste0("Sample: ", snakemake@wildcards[["sample"]], " | ASHLEYS probabilities"))

dev.off()

write.table(ashleys_data, file = snakemake@output[["well_table"]], sep = "\t", row.names = FALSE, quote = FALSE)
