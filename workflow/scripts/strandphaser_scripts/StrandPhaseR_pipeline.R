#!/usr/bin/Rscript
options(error = traceback)
args <- commandArgs(TRUE)

# add user defined path to load needed libraries
.libPaths(c(.libPaths(), args[6]))

# suppressPackageStartupMessages(library(StrandPhaseR))

# FIXME : tmp debuging local repo
library(devtools)

# load package w/o installing
# load_all("/g/korbel2/weber/Gits/StrandPhaseR")
strandphaser_path <- "/g/korbel2/weber/Gits/Strandphaser_clean/StrandPhaseR"
print(strandphaser_path)
load_all(strandphaser_path)

strandPhaseR(inputfolder = args[1], outputfolder = args[2], configfile = args[3], WCregions = args[4], positions = args[5], fillMissAllele = args[5])
