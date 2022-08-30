args <- commandArgs(TRUE)

# add user defined path to load needed libraries
# .libPaths(c(.libPaths(), args[6]))

# library(StrandPhaseR)
library(devtools)

# source("/g/korbel2/weber/Gits/StrandPhaseR/R/StrandPhase.R")
load_all("/g/korbel2/weber/Gits/StrandPhaseR/")


strandPhaseR(
    inputfolder = "/g/korbel/tsapalou/anaconda3/mosaicatcher-pipeline/POOLING_POOL1/HGSVCxpool1/all",
    outputfolder = "./TEST_Celia_chr4",
    configfile = "/g/korbel/tsapalou/anaconda3/mosaicatcher-pipeline/OUTPUT_pool1/strandphaser/HGSVCxpool1/StrandPhaseR.chr4.config",
    WCregions = "/g/korbel/tsapalou/anaconda3/mosaicatcher-pipeline/OUTPUT_pool1/strandphaser/HGSVCxpool1/strandphaser_input.txt",
    positions = "/g/korbel/tsapalou/anaconda3/mosaicatcher-pipeline/OUTPUT_pool1/snv_calls/HGSVCxpool1/chr4.vcf",
)
