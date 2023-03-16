# Whoeps, 07th Jan 2021
# Making a large overview over the results from the arbigent folder.
# I'm giving myself 2h to make this nice today.

# Input: callmatrix from clean_genotype.R
# Input: a csv from david from which to extract samplenames
# Output: a matrix with added entries:
#   Filter - Pass, NoReadsPass, MendelFail, FalsePositive, lowconf, AlwaysComplex
#   Mapability - percent?
#   nhom, nhet, nref, nnoreads, ncomplex
#   mendel 1/0

# Load libraries
library(stringr)
library(dplyr)
library(pheatmap)
library(matrixStats)
library(reshape2)
library(optparse)
source("workflow/scripts/arbigent/postprocess_helpers.R")


# INPUT INSTRUCTIONS
option_list <- list(
  make_option(c("-i", "--table"),
    type = "character", default = NULL,
    help = "res.csv", metavar = "character"
  ),
  make_option(c("-n", "--normal_names"),
    type = "character", default = T,
    help = "Should samplenames be converted from GM to NA?", metavar = "character"
  ),
  make_option(c("-o", "--outfile"),
    type = "character", default = NULL,
    help = "Outfile: verdicted table", metavar = "character"
  )
)

# Parse input
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
callmatrix_link <- opt$table
normal_names <- opt$normal_names
outfile <- opt$outfile

# callmatrix_link = '~/s/g/korbel2/StrandSeq/Test_WH/pipeline_7may/pipeline/regenotyper_allsamples_bulk/arbigent_results/res.csv'
# normal_names = T

cm <- read.table(callmatrix_link, header = 1, sep = "\t", stringsAsFactors = F)


### GO ###

# some inventory. Which samples do we have here? And therefore how many 'other' cols?
print(colnames(cm))
print(tail(colnames(cm), 1))

# samples <- colnames(cm)[grep("^[HMNG].*", colnames(cm))]
samples <- tail(colnames(cm), 1)
print(samples)
# stop()
n_samples <- length(samples)
n_other_cols <- dim(cm)[2] - n_samples

# Rename samples if wanted
if (as.numeric(normal_names)) {
  colnames(cm)[(n_other_cols + 1):dim(cm)[2]] <-
    str_replace(colnames(cm)[(n_other_cols + 1):dim(cm)[2]], "GM", "NA")
  # samples <- colnames(cm)[grep("^[HMNG].*", colnames(cm))]
  samples <- tail(colnames(cm), 1)
}

print(samples)
# Factor char stuff
cm[] <- lapply(cm, as.character)

# stratify entries with 0 valid bins.
cm[cm$valid_bins == 0, samples] <- "noreads"


# Count hom, het, ref, noreads and complex
cm <- count_homhetrefetc(cm, n_samples)

# Calc mapability

cm$valid_bins <- as.numeric(cm$valid_bins)

# Mendel
cm <- add_mendelfails(cm)

# Filter
cm <- apply_filter_new(cm, samples)

# Clean
cm[, c("mendel1", "mendel2", "mendel3")] <- NULL
cm[] <- lapply(cm, as.character)

# Rename inv_dup genotypes
cm <- make_invdups_human_readable(cm, samples)
print(cm)
# stop()
# Sort columns
cols <- c(colnames(cm)[1:n_other_cols], "verdict", "nref", "nhet", "nhom", "ninvdup", "ncomplex", samples)
cm_return <- cm[, cols]

# Save
write.table(cm_return, file = outfile, col.names = T, row.names = F, sep = "\t", quote = F)
