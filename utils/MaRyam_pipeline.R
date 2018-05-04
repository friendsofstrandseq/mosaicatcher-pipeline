#!/usr/bin/env Rscript

print(.libPaths())
sessionInfo()


args = commandArgs(trailingOnly=TRUE)
.libPaths( c( args[9],.libPaths()) )
suppressPackageStartupMessages(library(MaRyam))


args = as.data.frame(strsplit(args, split = "="), stringsAsFactors = F)

binRCfile = args[2,match("binRCfile", as.character(args[1,]))]
BRfile = args[2,match("BRfile", as.character(args[1,]))]
infoFile = args[2,match("infoFile", as.character(args[1,]))]
stateFile = args[2,match("stateFile", as.character(args[1,]))]
outputDir = args[2,match("outputDir", as.character(args[1,]))]
bin.size = as.numeric(args[2,match("bin.size", as.character(args[1,]))])
K = as.numeric(args[2,match("K", as.character(args[1,]))])
maximumCN = as.numeric(args[2,match("maximumCN", as.character(args[1,]))])

haplotypeMode=F
if (any(as.character(args[1,])=="haplotypeMode")){haplotypeMode = T}

print(paste("binRCfile =", binRCfile))
print(paste("BRfile =", BRfile))
print(paste("infoFile =", infoFile))
print(paste("stateFile =", stateFile))
print(paste("outputDir =", outputDir))
print(paste("bin.size =", bin.size))
print(paste("K =", K))
print(paste("maximumCN =", maximumCN))

p = data.table::fread(infoFile)$nb_p[1]

l <- changeRCformat(binRCfile, outputDir, p = p)
cellNames <- l$cellNames
initial.binRC <- l$binRC
r = matrix(rep(l$r, K), ncol = length(cellNames), byrow = T)
# report a warning if some cells are totally removed (because of having SCEs in all chrs)
if (length(r) != K*length(cellNames))
{
  message("Warning! The dimension of dispersion parameters and the number
          of cells and chromosomes don't match.
           Some cells may have been removed completely!")
}
f <- factor(initial.binRC$chromosome, levels=unique(initial.binRC$chromosome))
binRC <- split(initial.binRC, f)

cellTypes = changeCellTypesFormat(stateFile, cellNames)
#NBparams = changeNBparamsFormat(infoFile, K, cellNames)
#p = NBparams[[1]]
#r = NBparams[[2]]


segmentsCounts = getSegReadCounts(binRC, BRfile, K, bin.size)

SVcalling_wrapperFunc(bin.size, K, maximumCN, segmentsCounts, r, p, cellTypes, outputDir, haplotypeMode = haplotypeMode)

