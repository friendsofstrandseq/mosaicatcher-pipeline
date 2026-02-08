library(platetools)
library(ggplot2)
library(viridis)
library(dplyr)

## collect ASHLEYS prediction and count files
# ashleys_data <- read.table(file = "/scratch/tweber/DATA/TMP/labels384.tsv", sep = "\t", header = TRUE)
ashleys_data <- read.table(file = snakemake@input[["labels"]], sep = "\t", header = TRUE)
num_cells <- nrow(ashleys_data)
ashleys_data <- dplyr::arrange(ashleys_data, cell)
colnames(ashleys_data)[1] <- "ashleys_id"

# Extract well numbers from cell IDs
# Cell ID format: SampleNameiT[Barcode][WellNumber].sort.mdup.bam
# Example: LanexLorenzoNPCDMSO18holdSortUVLEDiTRU3C01.sort.mdup.bam -> 01
# Need to extract digits that come after a letter code and before .sort.mdup.bam

extract_well_number <- function(cell_id) {
    # Remove the file extension
    cell_clean <- sub("\\.sort\\.mdup\\.bam$", "", cell_id)

    # Extract last 2 digits from the cleaned cell ID
    matches <- regmatches(cell_clean, regexpr("\\d{2}$", cell_clean))
    if (length(matches) > 0 && nchar(matches[1]) > 0) {
        return(as.numeric(matches[1]))
    } else {
        return(NA_real_)
    }
}

# Determine plate type based on number of cells
if (num_cells > 0) {
    well_numbers <- vapply(ashleys_data$cell, extract_well_number, numeric(1))
    max_well <- max(well_numbers, na.rm = TRUE)

    # Estimate plate type from max well number
    if (max_well <= 96) {
        plate_type <- 96
    } else {
        plate_type <- 384
    }
} else {
    plate_type <- 96  # default
}

# Convert well number to well position (A1, A2, ..., H12 for 96-well, etc.)
Well_position <- character(length = num_cells)

for (idx in 1:num_cells) {
    well_num <- extract_well_number(ashleys_data$cell[idx])

    if (is.na(well_num)) {
        Well_position[idx] <- "Unknown"
    } else if (plate_type == 96) {
        # 96-well: 8 rows (A-H), 12 columns (1-12)
        row_idx <- ((well_num - 1) %% 8) + 1
        col_idx <- ((well_num - 1) %/% 8) + 1
        Well_position[idx] <- paste0(LETTERS[row_idx], col_idx)
    } else if (plate_type == 384) {
        # 384-well: 16 rows (A-P), 24 columns (1-24)
        row_idx <- ((well_num - 1) %% 16) + 1
        col_idx <- ((well_num - 1) %/% 16) + 1
        Well_position[idx] <- paste0(LETTERS[row_idx], col_idx)
    } else {
        Well_position[idx] <- "Unknown"
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
