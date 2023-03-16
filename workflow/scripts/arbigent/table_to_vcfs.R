# Whoeps, 10th Feb 2021.
# This is the main script for step three of the regenotyper Snakemake.
# Take an all.txt and turn it into a series of vcfs.
# Filtering and testing will, i think, be done by another file.

library(ggplot2)
library(reshape)
library(dplyr)
library(tibble)
library(optparse)

source("workflow/scripts/arbigent/clean_genotype_defs.R")
source("workflow/scripts/arbigent/vcf_defs.R")

# INPUT INSTRUCTIONS
option_list <- list(
  make_option(c("-a", "--alltxt"),
    type = "character", default = NULL,
    help = "ArbiGent output (traditionally 'all.txt') to be turned into vcf", metavar = "character"
  ),
  make_option(c("-m", "--msc"),
    type = "character", default = NULL,
    help = "count debug file from mosaicatcher main run", metavar = "character"
  ),
  make_option(c("-o", "--outdir"),
    type = "character", default = "./outputcorr/",
    help = "Outputdir", metavar = "character"
  ),
  make_option(c("-c", "--copynumberann"),
    type = "logical", default = F,
    help = "Use the copynumberannotation?", metavar = "character"
  )
)
# Parse input
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
alltxt_file <- opt$alltxt
msc_file <- opt$msc
outdir <- opt$outdir
use_cntrack <- opt$c
# alltxt_file = "/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/U32_freezemerge/sv_probabilities/all.txt"
# msc_file = "/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/U32_freezemerge/msc.debug"
# outdir = "/home/hoeps/Desktop/arbitrash"


### PARAMETERS ###

# Cutoff for 'lowconf'
save <- T
cutoff <- 3
# Second criterion: plus 5
# Complex vs simple: we want complex LLHs to be double the ones
# of simple, and at the same time at least higher by magnitude 5.
bias_add_factor <- 5


### FUNCTIONS ###
load_tab <- function(alltxt_file_f) {
  # Tiny function to load the input file
  tab <- read.table(alltxt_file_f, header = T, stringsAsFactors = F)
  tab <- tab %>% mutate(ID = paste0(chrom, "-", start + 1, "-INV-", (end - start) + 1))

  return(tab)
}

simplify_countmatrix_idup <- function(countm) {
  # Strongly simplify things

  countm[countm == "noreads"] <- "./."
  countm[countm == "1101"] <- "./1"
  countm[countm == "0100"] <- "1|."
  countm[countm == "2101"] <- "./1"
  countm[countm == "1110"] <- "./0"
  countm[countm == "0010"] <- ".|0"
  countm[countm == "0103"] <- "1|1"
  countm[countm == "1201"] <- "./1"
  countm[countm == "0001"] <- "./1"
  countm[countm == "1030"] <- "0|0"
  countm[countm == "0301"] <- "1|1"
  countm[countm == "3010"] <- "0|0"
  countm[countm == "0120"] <- "1|0"
  countm[countm == "3001"] <- "./1"
  countm[countm == "2020"] <- "0|0"
  countm[countm == "0202"] <- "1|1"
  countm[countm == "1000"] <- "0|."

  countm[countm == "1111"] <- "idup_1|1"
  countm[countm == "2200"] <- "idup_1|1"
  countm[countm == "0022"] <- "idup_1|1"

  countm[countm == "1011"] <- "idup_0|1"
  countm[countm == "1110"] <- "idup_1|0"

  simple_calls <- c(
    "0|0", "0|1", "1|0", "1|1",
    "0/0", "0/1", "1/0", "1/1",
    "./1", "1/.", "0/.", "./0",
    ".|1", "1|.", "0|.", ".|0"
  )

  # for idups, we want to make presence-absence variation analysis
  for (row in seq(1:dim(countm)[1])) {
    for (col in seq(1:dim(countm)[2])) {
      if ((countm[row, col] %in% simple_calls)) {
        countm[row, col] <- "0|0"
      }
    }
  }
  print("hi")
  good_calls <- c(
    "0|0", "0|1", "1|0", "1|1",
    "0/0", "0/1", "1/0", "1/1",
    "./1", "1/.", "0/.", "./0",
    ".|1", "1|.", "0|.", ".|0",
    "idup_0|1", "idup_1|0", "idup_1|1"
  )

  # Ugly loop because i dont know how to do it better
  for (row in seq(1:dim(countm)[1])) {
    for (col in seq(1:dim(countm)[2])) {
      if (!(countm[row, col] %in% good_calls)) {
        countm[row, col] <- "./."
      }
    }
  }
  head(countm)
  return(countm)
}

simplify_countmatrix <- function(countm) {
  # Strongly simplify things

  countm[countm == "noreads"] <- "./."
  countm[countm == "1101"] <- "./1"
  countm[countm == "0100"] <- "1|."
  countm[countm == "2101"] <- "./1"
  countm[countm == "1110"] <- "./0"
  countm[countm == "0010"] <- ".|0"
  countm[countm == "0103"] <- "1|1"
  countm[countm == "1201"] <- "./1"
  countm[countm == "0001"] <- "./1"
  countm[countm == "1030"] <- "0|0"
  countm[countm == "0301"] <- "1|1"
  countm[countm == "3010"] <- "0|0"
  countm[countm == "0120"] <- "1|0"
  countm[countm == "3001"] <- "./1"
  countm[countm == "2020"] <- "0|0"
  countm[countm == "0202"] <- "1|1"
  countm[countm == "1000"] <- "0|."

  simple_calls <- c(
    "0|0", "0|1", "1|0", "1|1",
    "0/0", "0/1", "1/0", "1/1",
    "./1", "1/.", "0/.", "./0",
    ".|1", "1|.", "0|.", ".|0"
  )
  # Ugly loop because i dont know how to do it better
  for (row in seq(1:dim(countm)[1])) {
    for (col in seq(1:dim(countm)[2])) {
      if (!(countm[row, col] %in% simple_calls)) {
        countm[row, col] <- "./."
      }
    }
  }
  head(countm)
  return(countm)
}

load_and_prep_CN <- function(msc_file_f) {
  CN <- read.table(msc_file, header = 1, stringsAsFactors = F)
  CNmerge <- as.data.frame(lapply(CN[, c("chrom", "start", "end", "valid_bins")], as.character))

  CNmerge <- tibble::as_tibble(CNmerge)
  CNmerge <- CNmerge %>% mutate(
    chrom = as.character(chrom),
    start = as.numeric(as.character(start)),
    end = as.numeric(as.character(end)),
    valid_bins = as.numeric(as.character(valid_bins))
  )
  return(CNmerge)
}


################
### RUN CODE ###
################

# Load input file
tab <- load_tab(alltxt_file)
# Load CN file. It will be used to include 'valid bins' information, which is good to have in the output files.
CNmerge <- load_and_prep_CN(msc_file)

# TODO: WHAT IS THIS LINE DOING?
bias_factor <- tab$confidence_hard_over_second # This is the first criterion: double


######### GET 'REPORTED' GENOTYPES ##############
# tabp = [tab]le_[p]rocessed. This has added a 'GT' column that is the GT result
# of choice. And this GT of choice depends on: bias(add)factor, cutoff and wether
# or not we want 'lowconf' label included.

if (use_cntrack) {
  tab$pred_hard <- paste(tab$pred_hard, tab$illumina_CN, sep = ":")
  tab$pred_nobias <- paste(tab$pred_nobias, tab$illumina_CN, sep = ":")
}
# tabp = Complex calls allowed, lowconf label added
tabp <- add_gts_revisited_lowconf(tab, bias_factor, bias_add_factor, cutoff)
print(tabp)
# tabp2 =  Complex calls allowed, lowconf label nope, LLHs printed
tabp2 <- add_long_gts_revisited(tab, bias_factor, bias_add_factor, cutoff)
print(tabp2)

# tabp3 = Complex calls allowed, lowconf label nope
tabp3 <- add_gts_revisited(tab, bias_factor, bias_add_factor, cutoff)


# Merge valid bin information into tabp's
tabp <- full_join(tabp, CNmerge, by = c("chrom", "start", "end"))
tabp2 <- full_join(tabp2, CNmerge, by = c("chrom", "start", "end"))
tabp3 <- full_join(tabp3, CNmerge, by = c("chrom", "start", "end"))
print(tabp3)
print(CNmerge)

####################################################################################

# Cast table tabp3 into vcf-like matrix
callmatrix <- cast(unique(tabp3), chrom + start + end + ID + len + valid_bins ~ sample, value = "GT")

####################################################################################

# Name shortening, a bit of reformatting
cms <- callmatrix
print(colnames(cms))
samplenames <- colnames(cms)[7:length(colnames(cms))]
print(samplenames)

cms_gts <- cms[, samplenames, drop = F]
print(cms_gts)

# print(head(lapply(cms_gts, as.character)))
cms_gts[] <- lapply(cms_gts, as.character)
# Complex calls to simple ones
print(cms_gts)

countm <- simplify_countmatrix(cms_gts)
countm_idup <- simplify_countmatrix_idup(cms_gts)
# Bind description and GTs back together
cms_full <- cbind(cms[, 1:6], cms_gts)

# Bind simplified countmatrix to description
cms_simple <- cbind(cms[, 1:6], countm)
cms_simple_idup <- cbind(cms[, 1:6], countm_idup)
####################################################################################

# Tabp: this is the 'normal' one. With lowconf label
callmatrix <- cast(unique(tabp), chrom + start + end + ID + len + valid_bins ~ sample, value = "GT")
callmatrix <- callmatrix[, colSums(is.na(callmatrix)) < dim(callmatrix)[2]]

# Tabp2: this is the one with LLH info included
callmatrix_detail <- cast(tabp2, chrom + start + end + ID + len + valid_bins ~ sample, value = "GTL")
callmatrix_detail <- callmatrix_detail[, colSums(is.na(callmatrix_detail)) < dim(callmatrix_detail)[2]]
# callmatrix_detail = callmatrix_detail[ , colSums(is.na(callmatrix_detail)) == 0]

# Sidequest: find hom invs. A bit messy but we keep it for now.
callmatrix_hom_lab <- callmatrix
callmatrix_hom_lab$nhom <- rowSums(callmatrix_hom_lab == "1|1")
cm_hom <- callmatrix_hom_lab[callmatrix_hom_lab$nhom > (dim(callmatrix_hom_lab)[2] - 5) * 0.8, ]
hom_ids <- cm_hom$ID
cm_detail_hom <- callmatrix_detail[callmatrix_detail$ID %in% hom_ids, ]

################## MAKE VCFS
# The detailed one from tabp2, including LLH info.
vcf <- vcfify_callmatrix_detail(callmatrix_detail)
# Like above, but filtered for misos
vcf_miso <- vcfify_callmatrix_detail(cm_detail_hom)
# And this one is based on simplified tabp3.
vcf_limix <- vcfify_callmatrix_simple_for_limix(cms_simple)
# And here we save the one that contains idups
vcf_limix_plus_idups <- vcfify_callmatrix_simple_for_limix(cms_simple_idup)
################## Make at least one plot TODO
print(head(cms_simple_idup))
# ggplot(tabp3) + geom_point(aes(x=log1p(confidence_hard_over_second), y=log1p(confidence_nobias_over_hard), col=GT))


################# SAVE ALL THESE DIFFERENT THINGS.
if (save == T) {
  # Prep directory
  dir.create(outdir)

  # Paths, paths, paths.
  callmatrix_file <- "res.csv"
  callmatrix_file_detail <- "res_detail.csv"
  vcffile_all <- "res_all.vcf"
  vcffile_miso <- "res_miso.vcf"
  vcffile_limix <- "res_verysimple.vcf"
  vcffile_limix_plus_idups <- "res_verysimple_idups.vcf"

  # Save simple callmatrix
  write.table(callmatrix, file = file.path(outdir, callmatrix_file), quote = F, col.names = T, row.names = F, sep = "\t")
  # Same detailed callmatrix
  write.table(callmatrix_detail, file = file.path(outdir, callmatrix_file_detail), quote = F, col.names = T, row.names = F, sep = "\t")
  # Save vcf_all
  outvcffile <- file.path(outdir, vcffile_all)
  writeLines(vcf[[1]], file(outvcffile))
  write.table(vcf[[2]], file = outvcffile, quote = F, col.names = F, row.names = F, sep = "\t", append = T)
  system(paste0("cat ", outvcffile, ' | awk \'$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}\' > ', outvcffile, "_sorted"))

  # Save vcf_miso
  outvcffile_miso <- file.path(outdir, vcffile_miso)
  writeLines(vcf_miso[[1]], file(outvcffile_miso))
  write.table(vcf_miso[[2]], file = outvcffile_miso, quote = F, col.names = F, row.names = F, sep = "\t", append = T)
  system(paste0("cat ", outvcffile_miso, ' | awk \'$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}\' > ', outvcffile_miso, "_sorted"))

  # Save vcf_limix
  outvcffile_limix <- file.path(outdir, vcffile_limix)
  writeLines(vcf_limix[[1]], file(outvcffile_limix))
  write.table(vcf_limix[[2]], file = outvcffile_limix, quote = F, col.names = F, row.names = F, sep = "\t", append = T)
  system(paste0("cat ", outvcffile_limix, ' | awk \'$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}\' > ', outvcffile_limix, "_sorted"))

  # Save vcf_limix idups
  outvcffile_limix_idups <- file.path(outdir, vcffile_limix_plus_idups)
  writeLines(vcf_limix_plus_idups[[1]], file(outvcffile_limix_idups))
  write.table(vcf_limix_plus_idups[[2]], file = outvcffile_limix_idups, quote = F, col.names = F, row.names = F, sep = "\t", append = T)
  system(paste0("cat ", outvcffile_limix_idups, ' | awk \'$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}\' > ', outvcffile_limix_idups, "_sorted"))
}
# Save only complex stuff Complex exploration
# tabcomp = tab[tab$GTs == 'complex',]
# indiv_invs <- unique( tabcomp[ , 1:3 ] )
# library(dplyr)
# aa = left_join(indiv_invs,callmatrix_detail)
# bb = aa[
#  with(aa,
