args <- commandArgs(TRUE)

# add user defined path to load needed libraries
# .libPaths(c(.libPaths(), args[6]))

# library(StrandPhaseR)
library(devtools)

# source("/g/korbel2/weber/Gits/StrandPhaseR/R/StrandPhase.R")
load_all("/g/korbel2/weber/Gits/StrandPhaseR/")


strandPhaseR(
    inputfolder = "/g/korbel2/weber/MosaiCatcher_files/HGSVC_WH/GM12329/all",
    outputfolder = "/g/korbel2/weber/MosaiCatcher_output/HGSVC_WH/strandphaser/GM12329/StrandPhaseR_analysis.chr21",
    configfile = "/g/korbel2/weber/MosaiCatcher_output/HGSVC_WH/strandphaser/GM12329/StrandPhaseR.chr21.config",
    WCregions = "/g/korbel2/weber/MosaiCatcher_output/HGSVC_WH/strandphaser/GM12329/strandphaser_input.txt",
    positions = "/g/korbel2/weber/MosaiCatcher_output/HGSVC_WH/snv_genotyping/GM12329/chr21.vcf",
)