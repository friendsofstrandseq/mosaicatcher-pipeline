#!/usr/bin/Rscript
options(error = traceback)
args <- commandArgs(TRUE)

# add user defined path to load needed libraries
.libPaths(c(.libPaths(), args[6]))

suppressPackageStartupMessages(library(StrandPhaseR))

# FIXME : tmp debuging local repo
# library(devtools)

# load package w/o installing
# load_all("/g/korbel2/weber/Gits/StrandPhaseR")
# strandphaser_path <- "/g/korbel2/weber/Gits/Strandphaser_clean/StrandPhaseR"
# print(strandphaser_path)
# load_all(strandphaser_path)

# --- Robustness wrapper (does NOT change results for well-phased chromosomes) ---
# StrandPhaseR's fillGapsWithVCF crashes on degenerate input (a chromosome with too
# few informative SNPs collapses to a length-1 GRanges). That crash kills the whole
# job, and because combine_strandphaser_output needs every chromosome, a single
# sparse chromosome blocks the entire sample. We can't patch the packaged StrandPhaseR,
# so we catch the failure here and emit the two expected output files empty
# (header-only) so the chromosome is simply "not phased" instead of fatal.
#
# Consistency: write_placeholder_outputs() only writes files that are MISSING or empty.
# On a successful chromosome StrandPhaseR has already written non-empty outputs, so this
# is a no-op there and those results are untouched.

# Derive sample + chrom from the paths (used only to build placeholder outputs).
# args[1] = <data_location>/<sample>/selected ; args[2] = <...>/StrandPhaseR_analysis.<chrom>
sample <- basename(dirname(args[1]))
chrom <- sub("^.*StrandPhaseR_analysis\\.", "", args[2])

write_placeholder_outputs <- function() {
    phased_dir <- file.path(args[2], "Phased")
    vcf_dir <- file.path(args[2], "VCFfiles")
    dir.create(phased_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(vcf_dir, recursive = TRUE, showWarnings = FALSE)

    phased_haps <- file.path(phased_dir, "phased_haps.txt")
    if (!file.exists(phased_haps) || file.info(phased_haps)$size == 0) {
        writeLines(
            "sample\tcell\tchrom\tstart\tend\tclass\thap1.cis.simil\thap1.trans.simil\thap2.cis.simil\thap2.trans.simil",
            phased_haps
        )
    }

    vcf_file <- file.path(vcf_dir, paste0(chrom, "_phased.vcf"))
    if (!file.exists(vcf_file) || file.info(vcf_file)$size == 0) {
        writeLines(
            c(
                "##fileformat=VCFv4.2",
                "##source=StrandPhase_algorithm",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                paste0("##contig=<ID=", chrom, ">"),
                paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", sample)
            ),
            vcf_file
        )
    }
}

tryCatch(
    strandPhaseR(inputfolder = args[1], outputfolder = args[2], configfile = args[3], WCregions = args[4], positions = args[5], fillMissAllele = args[5]),
    error = function(e) {
        message(
            "[WARNING] StrandPhaseR failed for ", sample, " ", chrom,
            " (likely too few informative SNPs to phase): ", conditionMessage(e)
        )
        message("[WARNING] Writing empty placeholder outputs so this chromosome is skipped without blocking the sample.")
    }
)

# Guarantee the two declared output files exist (no-op when StrandPhaseR succeeded).
write_placeholder_outputs()
