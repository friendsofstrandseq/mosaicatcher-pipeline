options(error = traceback)
args <- commandArgs(TRUE)

# add user defined path to load needed libraries
.libPaths(c(.libPaths(), args[6]))

# library(StrandPhaseR)
library(devtools)

# source("/g/korbel2/weber/Gits/StrandPhaseR/R/StrandPhase.R")
load_all("/g/korbel2/weber/Gits/StrandPhaseR/")
print("/g/korbel2/weber/Gits/StrandPhaseR/")

# suppressPackageStartupMessages(library(StrandPhaseR))


strandPhaseR(
    inputfolder = args[1],
    outputfolder = args[2],
    configfile = args[3],
    WCregions = args[4],
    positions = args[5],
)
