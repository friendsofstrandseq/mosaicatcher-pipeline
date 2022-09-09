options(error = traceback)
args <- commandArgs(TRUE)

# add user defined path to load needed libraries
.libPaths(c(.libPaths(), args[6]))

# library(StrandPhaseR)
# library(devtools)

# source("/g/korbel2/weber/Gits/StrandPhaseR/R/StrandPhase.R")
# load_all("/g/korbel2/weber/Gits/StrandPhaseR/")


suppressPackageStartupMessages(library(StrandPhaseR))


strandPhaseR(
    inputfolder = "/g/korbel2/weber/MosaiCatcher_files/POOLING/PSEUDOPOOL/pseudopool/all",
    outputfolder = "./TEST_SPHR",
    configfile = "/g/korbel/weber/MosaiCatcher_output/POOLING/PSEUDOPOOL/strandphaser/pseudopool/StrandPhaseR.chr17.config",
    WCregions = "/g/korbel/weber/MosaiCatcher_output/POOLING/PSEUDOPOOL/strandphaser/pseudopool/strandphaser_input.txt",
    positions = "/g/korbel2/weber/workspace/mosaicatcher-update/test_chr17.vcf",
)
