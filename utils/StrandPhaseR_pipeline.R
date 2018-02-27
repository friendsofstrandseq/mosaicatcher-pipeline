#!/usr/bin/Rscript

args=commandArgs(TRUE)

#add user defined path to load needed libraries
.libPaths( c( .libPaths(), args[6]) )

suppressPackageStartupMessages(library(StrandPhaseR))

strandPhaseR(inputfolder=args[1], outputfolder=args[2], configfile = args[3],  WCregions = args[4] , positions=args[5])
