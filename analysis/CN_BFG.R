brs <- c(33998805, 34185186, 34211372, 35497754, 37452507, 39578189)
breakpoints <- data.table(k = length(unique(round(brs/100000))), chrom = "chr10", bps = unique(round(brs/100000)))
brfile <- "/home/maryam/research/hackathons/C7-7Jun/segmentation2/C7/manual-regions.txt"
write.table(breakpoints, file = brfile, quote = F, sep = "\t", row.names = F)

dir <- "/home/maryam/research/hackathons/C7-7Jun"
countsFile <- file.path(dir, "100000_fixed.txt.gz")
brFile <- file.path(dir, "segmentation2/C7/manual-regions.txt")
infoFile <- file.path(dir, "100000_fixed.info")
stateFile <- file.path(dir, "final.txt")

# counts[chrom=="chr10" & start > 33998805 & end < 34185186]
counts <- fread(paste("zcat", countsFile))
segs <- fread(brFile)
info <- fread(infoFile)
strand <- fread(stateFile)

counts[class=="None", class:="WW"]
probs <- mosaiClassifierPrepare(counts, info, strand, segs)
probs[, nb_r:=expected*nb_p/(1-nb_p)]

maximumCN <- 200

probs <- probs[, cbind(.SD, CN=1:maximumCN), 
               by=.(sample, cell, chrom, start, end)]


CNprobs <- probs[, .(cov=sum(W+C), 
                     expected=sum(expected), 
                     nb_p=nb_p[1], 
                     nb_r=sum(nb_r)), 
                 by=.(sample, chrom, start, end)]
# first observation
(CNprobs$cov / CNprobs$expected) * 2

CNprobs <- CNprobs[, cbind(.SD, CN=1:maximumCN), 
               by=.(sample, chrom, start, end)]


# computing NB probs and calling most probable CN in probs table per segment per cell
probs[, CN_ll:=dnbinom(x=C+W, size=CN*(nb_r/2), prob=nb_p)]
cellCNcalls <- probs[, .SD[which.max(CN_ll)], 
        by=.(sample, cell, chrom, start, end)]

# computing NB probs and calling CN in CNprobs table per segment
CNprobs[, CN_ll:=dnbinom(x=cov, size=CN*(nb_r/2), prob=nb_p)]
CNcalls <- CNprobs[, .SD[which.max(CN_ll)], 
        by=.(sample, chrom, start, end)]


# removing extra columns for cleaning the data for plotting phylogenetic tree
cellCNcalls <- cellCNcalls[start!=0]
cellCNcalls[, `:=`(nb_p=NULL, class=NULL, C=NULL, W=NULL,expected=NULL, scalar=NULL, nb_r=NULL, CN_ll=NULL)]
write.table(cellCNcalls, file.path(dir, "BFB_cell_CNs.table"))
