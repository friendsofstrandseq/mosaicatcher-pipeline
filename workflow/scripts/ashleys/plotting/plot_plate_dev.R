library(platetools)
library(ggplot2)
library(viridis)
library(dplyr)
library(stringr)


# df <- data.frame(vals = rnorm(384),
#                  well = num_to_well(1:384, plate = 384))

# print(df)
# stop()

args <- commandArgs(trailingOnly = T)

prefix = args[1]

## collect ASHLEYS prediction and count files
ashleys_data <- read.table(file =  paste0(prefix, "/cell_selection/labels.tsv"), sep = "\t", header = TRUE)
# ashleys_data <- read.table(file = snakemake@input[["labels"]], sep = "\t", header = TRUE)
plate_type <- nrow((ashleys_data))
ashleys_data <- dplyr::arrange(ashleys_data, cell)
colnames(ashleys_data)[1] <- "ashleys_id"

corr_table_path = "workflow/data/plotting/384_1A3C5E7G_correspondance_table.tsv"
corr_table <- read.table(corr_table_path, header = TRUE, sep = "\t")
sample <- basename(prefix)


# Apply regex and extract groups
ashleys_data <- ashleys_data %>%
  mutate(
    index = str_extract(ashleys_id, "(iTRU|PE20)[A-Za-z0-9]{4,5}")
  )


# View the result

# Well_position <- character()

# if (plate_type == 96) {
#     for (i in 1:12)
#     {
#         for (j in 1:8)
#         {
#             tmp <- paste0(LETTERS[j], i)
#             Well_position <- c(Well_position, tmp)
#         }
#     }
# } else if (plate_type == 384) {
#     for (i in 1:24)
#     {
#         for (j in 1:16)
#         {
#             tmp <- paste0(LETTERS[j], i)
#             Well_position <- c(Well_position, tmp)
#         }
#     }
# }


print(corr_table)
print(ashleys_data)
ashleys_data <- merge(ashleys_data, corr_table, by.x = "index", by.y = "index", all.x = TRUE)

write.table(ashleys_data, file = paste0(prefix, "/plots/plate/ashleys_well_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


ashleys_data <- ashleys_data %>%
  mutate(
    Well_row = str_extract(Well_position, "[A-Za-z]+"),
    Well_col = as.integer(str_extract(Well_position, "\\d+"))
  ) %>%
  arrange(Well_row, Well_col)


print(ashleys_data)


pdf(paste0(prefix, "/plots/plate/ashleys_plate_predictions.pdf"))
# pdf(snakemake@output[["predictions"]])


raw_map(
    data = ashleys_data$prediction,
    well = ashleys_data$Well_position,
    plate = plate_type
) +
    scale_fill_distiller(type = "div", palette = "RdYlGn", direction = 1) +
    # ggtitle(paste0("Sample: TEST | ASHLEYS binary predictions (cutoff=0.5)"))
    ggtitle(paste0("Sample: ", sample, " | ASHLEYS binary predictions (cutoff=", 0.5, ")"))

dev.off()

# pdf("TEST_ashleys_plate_probabilities.pdf")
pdf(paste0(prefix, "/plots/plate/ashleys_plate_probabilities.pdf"))

# pdf(snakemake@output[["probabilities"]])

raw_map(
    data = ashleys_data$probability,
    well = ashleys_data$Well_position,
    plate = plate_type
) +
    scale_fill_distiller(type = "div", palette = "RdYlGn", direction = 1) +
    # ggtitle(paste0("Sample: ", "TEST", " | ASHLEYS probabilities"))
    ggtitle(paste0("Sample: ", sample, " | ASHLEYS probabilities"))

dev.off()