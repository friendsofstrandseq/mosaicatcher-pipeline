args <- commandArgs(TRUE)

# add user defined path to load needed libraries
# .libPaths(c(.libPaths(), args[6]))

# library(StrandPhaseR)
library(devtools)

# source("/g/korbel2/weber/Gits/StrandPhaseR/R/StrandPhase.R")
load_all("/g/korbel2/weber/Gits/StrandPhaseR/")


strandPhaseR(
    inputfolder = "/g/korbel2/weber/MosaiCatcher_files/LCL_TALL/H2NCTAFX2_GM20509B_20s000579-1-1/selected",
    outputfolder = "/g/korbel2/weber/MosaiCatcher_output/Mosaicatcher_output_singularity_LCL-TALL/strandphaser/H2NCTAFX2_GM20509B_20s000579-1-1/StrandPhaseR_analysis.chrX",
    configfile = "/g/korbel2/weber/MosaiCatcher_output/Mosaicatcher_output_singularity_LCL-TALL/strandphaser/H2NCTAFX2_GM20509B_20s000579-1-1/StrandPhaseR.chrX.config",
    WCregions = "/g/korbel2/weber/MosaiCatcher_output/Mosaicatcher_output_singularity_LCL-TALL/strandphaser/H2NCTAFX2_GM20509B_20s000579-1-1/strandphaser_input.txt",
    positions = "/g/korbel2/weber/MosaiCatcher_output/Mosaicatcher_output_singularity_LCL-TALL/snv_genotyping/H2NCTAFX2_GM20509B_20s000579-1-1/chrX.vcf",
)