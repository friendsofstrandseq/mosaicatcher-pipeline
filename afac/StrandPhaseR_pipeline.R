options(error = traceback)
args <- commandArgs(TRUE)

# add user defined path to load needed libraries
.libPaths(c(.libPaths(), args[6]))

# library(StrandPhaseR)
library(devtools)

# source("/g/korbel2/weber/Gits/StrandPhaseR/R/StrandPhase.R")
load_all("/g/korbel2/weber/Gits/StrandPhaseR/")


# suppressPackageStartupMessages(library(StrandPhaseR))


strandPhaseR(
    inputfolder = "/g/korbel2/weber/MosaiCatcher_files/POOLING/POOLING_POOL2_HG00186/HGSVCxpool2/all",
    outputfolder = "./TEST_SPHR",
    configfile = "/g/korbel2/weber/MosaiCatcher_output/POOLING/POOLING_POOL2_HG00186/strandphaser/HGSVCxpool2/StrandPhaseR.chr22.config",
    WCregions = "/g/korbel2/weber/MosaiCatcher_output/POOLING/POOLING_POOL2_HG00186/strandphaser/HGSVCxpool2/strandphaser_input.txt",
    positions = "/g/korbel2/weber/MosaiCatcher_output/POOLING/POOLING_POOL2_HG00186/snv_calls/HGSVCxpool2/chr22.vcf",
)
