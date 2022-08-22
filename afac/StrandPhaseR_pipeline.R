args <- commandArgs(TRUE)

# add user defined path to load needed libraries
# .libPaths(c(.libPaths(), args[6]))

# library(StrandPhaseR)
library(devtools)

# source("/g/korbel2/weber/Gits/StrandPhaseR/R/StrandPhase.R")
load_all("/g/korbel2/weber/Gits/StrandPhaseR/")


strandPhaseR(
    inputfolder = "/g/korbel2/weber/MosaiCatcher_files/POOLING/POOLING_POOL1/HGSVCxpool1/all",
    outputfolder = "/g/korbel2/weber/MosaiCatcher_output/POOLING/POOLING1_190822_100KB/strandphaser/HGSVCxpool1/StrandPhaseR_analysis.22",
    configfile = "/g/korbel2/weber/MosaiCatcher_output/POOLING/POOLING1_190822_100KB/strandphaser/HGSVCxpool1/StrandPhaseR.chr22.config",
    WCregions = "/g/korbel2/weber/MosaiCatcher_output/POOLING/POOLING1_190822_100KB/strandphaser/HGSVCxpool1/strandphaser_input.txt",
    positions = "/g/korbel2/weber/MosaiCatcher_output/POOLING/POOLING1_190822_100KB/snv_calls/HGSVCxpool1/chr22.vcf",
)
