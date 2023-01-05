
## from https://gist.github.com/Swarchal/b938933ae9ded94b3c14d6485b27cf69

## Requirements and library loading

# install.packages(c('platetools', 'ggplot2', 'viridis', 'dplyr'))

library(platetools)
library(ggplot2)
library(viridis)
library(dplyr)


## collect ASHLEYS prediction and count files
# prediction_file = list.files('.', 'labels_raw.tsv')

## read in ASHLEYS data
ashleys_data = read.table(file = "labels_test.tsv" , sep = '\t', header = TRUE)
# ashleys_data = read.table(file = snakemake@input[["labels"]] , sep = '\t', header = TRUE)
ashleys_data = dplyr::arrange(ashleys_data,cell)
## rename cell column because same as pipeline data
colnames(ashleys_data)[1] = 'ashleys_id'


## Create Well_position vector for plotting
Well_position = character()

for (i in 1:12)
{
    for (j in 1:8)
    {
        tmp = paste0(LETTERS[j], i)
        # print(tmp)
        Well_position = c(Well_position, tmp)
    }
}




### Plot 4: ASHLEYS prediction distribution


pdf('ASHLEYS_predictions.pdf')
# pdf(snakemake@output[["predictions"]])
ashleys_data$Well_position = Well_position
print(ashleys_data)

raw_map(data = ashleys_data$prediction,
       well = ashleys_data$Well_position, 
       plate = 96) +
    scale_fill_distiller(type = "div", palette = "RdYlGn", direction = 1)
    # ggtitle(paste0("Sample: ", snakemake@wildcards[["sample"]], " | ASHLEYS binary predictions (cutoff=", snakemake@config[["ashleys_threshold"]], ")"))

dev.off()

pdf('ASHLEYS_probabilities.pdf')
# pdf(snakemake@output[["probabilities"]])
ashleys_data$Well_position = Well_position
# print(ashleys_data)

raw_map(data = ashleys_data$probability,
       well = ashleys_data$Well_position, 
       plate = 96) +
    scale_fill_distiller(type = "div", palette = "RdYlGn", direction = 1)
    # ggtitle(paste0("Sample: ", snakemake@wildcards[["sample"]], " | ASHLEYS probabilities"))

dev.off()
