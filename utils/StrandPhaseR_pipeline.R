#!/usr/bin/Rscript

args=commandArgs(TRUE)

#add user defined path to load needed libraries
.libPaths( c( .libPaths(), args[6]) )

suppressPackageStartupMessages(library(StrandPhaseR))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))

strandPhaseR(inputfolder=args[1], outputfolder=args[2], configfile = args[3],  WCregions = args[4] , positions=args[5], numCPU=args[7])
