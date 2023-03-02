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
    inputfolder = "/g/korbel2/tsapalou/KAPA_POOL1/HGSVCpool1_15ulOP_96of384_KAPA/bam",
    outputfolder = "./TEST_KAPPA",
    configfile = "/g/korbel2/tsapalou/KAPA_POOL1/HGSVCpool1_15ulOP_96of384_KAPA/strandphaser/StrandPhaseR.chr5.config",
    WCregions = "/g/korbel2/tsapalou/KAPA_POOL1/HGSVCpool1_15ulOP_96of384_KAPA/strandphaser/strandphaser_input.txt",
    positions = "/g/korbel2/tsapalou/KAPA_POOL1/HGSVCpool1_15ulOP_96of384_KAPA/snv_calls/chr5.vcf",
)
