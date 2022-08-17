## load required libraries
library(GenomicAlignments)
library(ggplot2)
library(cowplot)
library(BiocParallel)

#' Print haplotagged read counts
#'
#' This function will take \code{list} of haplotagged bams files and will return \link{data.frame}
#' counts of reads per haplotype.
#'
#' @param sv.table A path to a table in bed format with regions to count haplotagged reads.
#' @param bam.path A path to the haplotagged bam files.
#' @author David Porubsky

getHaplotagTable <- function(sv.table = NULL, bam.path = NULL) {
  ## read the SV table
  SV.regions <- read.table(sv.table, header = TRUE, stringsAsFactors = FALSE)
  SV.regions.gr <- GRanges(seqnames = SV.regions$chrom, ranges = IRanges(start = SV.regions$start, end = SV.regions$end), cell = SV.regions$cell, class = SV.regions$class, sv_call_name = SV.regions$sv_call_name, sv_call_name_2nd = SV.regions$sv_call_name_2nd)
  ## list all bam files to count haplotagged reads in
  haplotag.bams <- list.files(path = bam.path, pattern = "\\.bam$", full.names = T)

  all.counts <- list()
  for (i in 1:length(haplotag.bams)) {
    bam <- haplotag.bams[i]
    filename <- basename(bam)
    message("Processing bamfile ", filename, " ...", appendLF = F)
    ptm <- proc.time()

    ## get regions for a give
    cell.id <- unlist(strsplit(filename, "\\."))[1]
    SV.regions.gr.perCell <- SV.regions.gr[SV.regions.gr$cell == cell.id]

    ## read in reads for selected regions
    fragments <- bamregion2GRanges(bamfile = bam, region = SV.regions.gr.perCell, pairedEndReads = T, min.mapq = 10, filterAltAlign = TRUE)

    fragments$HP[is.na(fragments$HP)] <- 0 # set missing haplotag to zero
    ## split reads per selected region
    hits <- findOverlaps(SV.regions.gr.perCell, fragments)
    fragments.per.region <- split(fragments[subjectHits(hits)], queryHits(hits))

    ## count haplotagged reads in selected regions
    counts <- lapply(fragments.per.region, getHapReadCount)
    counts.df <- do.call(rbind, counts)

    ## export final table of haplotagged read counts
    SV.regions.df <- as.data.frame(SV.regions.gr.perCell[unique(queryHits(hits))])[, c(1, 2, 3, 6, 7, 8, 9)]
    colnames(SV.regions.df) <- c("chrom", "start", "end", "cell", "class", "sv_call_name", "sv_call_name_2nd")
    cell.hap.counts <- cbind(SV.regions.df, counts.df)
    all.counts[[i]] <- cell.hap.counts

    time <- proc.time() - ptm
    message(" ", round(time[3], 2), "s")
  }

  final.table <- do.call(rbind, all.counts)
  file.base <- gsub(sv.table, pattern = "\\.txt|\\.bed|\\.csv|\\.tsv", replacement = "")
  file.destination <- paste0(file.base, "_haplotaggedCounts.txt")
  write.table(final.table, file = file.destination, quote = FALSE, row.names = FALSE)
  message("DONE!!!")

  return(final.table)
}

#' Print haplotagged read counts
#'
#' This function will take \code{list} of haplotagged bams files and will return \link{data.frame}
#' counts of reads per haplotype.
#'
#' @param bedFile A path to a table in bed format with regions to count haplotagged reads.
#' @param bam.file BAM file name
#' @param CPUs Number of CPUs to use for data processing.
#' @author David Porubsky

getHaplotagTable2 <- function(bedFile = NULL, bam.file = NULL, CPUs = 4, paired_end = TRUE, file.destination = NULL) {
  suppressPackageStartupMessages({
    requireNamespace("tools")
  })

  message("Creating haplotag table")
  message("BED file: ", bedFile)
  message("BAM file: ", bam.file)
  message("Output file: ", file.destination)

  ## read the SV table
  regions <- read.table(bedFile, header = FALSE, stringsAsFactors = FALSE)
  print(head(regions))
  print(dim(regions))
  regions.gr <- GRanges(seqnames = regions$V1, ranges = IRanges(start = regions$V2, end = regions$V3))
  print(regions.gr)

  filename <- basename(bam.file)
  cell.id <- unlist(strsplit(filename, "\\."))[1]
  message("Processing bamfile ", filename, " ...", appendLF = F)
  ptm <- proc.time()

  print(filename)
  print(cell.id)
  ## read in reads for selected regions

  fragments <- bamregion2GRanges(bamfile = bam.file, region = regions.gr, pairedEndReads = paired_end, min.mapq = 10, filterAltAlign = TRUE)



  fragments$HP[is.na(fragments$HP)] <- 0 # set missing haplotag to zero
  ## split reads per selected region
  hits <- findOverlaps(regions.gr, fragments)
  fragments.per.region <- split(fragments[subjectHits(hits)], queryHits(hits))
  # subset regions to the regions that have non-zero read count
  regions <- regions[unique(queryHits(hits)), ]
  ## count haplotagged reads in selected regions
  # use parallel execution with a given number of CPUs
  counts <- bplapply(fragments.per.region, getHapReadCount, BPPARAM = MulticoreParam(CPUs))
  counts.df <- do.call(rbind, counts)

  # cbind regions to counts.df
  counts.df <- cbind(cell = cell.id, regions, counts.df)
  # renaming the regions columns
  colnames(counts.df)[2:4] <- c("chrom", "start", "end")

  time <- proc.time() - ptm
  message(" ", round(time[3], 2), "s")

  if (!is.null(file.destination)) {
    write.table(counts.df, file = file.destination, quote = FALSE, row.names = FALSE)
  }
  message("DONE!!!")

  return(counts.df)
}

#' Count haplotype specific reads
#'
#' Get counts of haplotagged reads per strand stored \code{\link{GRanges}} object.
#'
#' @param gr A \code{\link{GRanges}} object.
#' @author David Porubsky

getHapReadCount <- function(gr) {
  ## get total watson and crick counts
  crick.reads <- gr[strand(gr) == "+"]
  watson.reads <- gr[strand(gr) == "-"]
  crick.count <- length(crick.reads)
  watson.count <- length(watson.reads)

  ## get read counts per haplotype and per directionality
  crick.H1 <- length(crick.reads[crick.reads$HP == 1])
  crick.H2 <- length(crick.reads[crick.reads$HP == 2])
  watson.H1 <- length(watson.reads[watson.reads$HP == 1])
  watson.H2 <- length(watson.reads[watson.reads$HP == 2])

  df.count <- data.frame(crick.count = crick.count, watson.count = watson.count, crick.H1 = crick.H1, crick.H2 = crick.H2, watson.H1 = watson.H1, watson.H2 = watson.H2)
  return(df.count)
}


#' Import BAM file into GRanges
#'
#' Import aligned reads from a BAM file into a \code{\link{GRanges}} object.
#'
#' @param file Bamfile with aligned reads.
#' @param bamindex Bam-index file with or without the .bai ending. If this file does not exist it will be created and a warning is issued.
#' @param region If only a subset of the genomic regions should be loaded.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param min.mapq Minimum mapping quality when importing from BAM files.
#' @importFrom Rsamtools indexBam scanBamHeader ScanBamParam scanBamFlag testPairedEndBam
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments first last
#' @author David Porubsky
#' @export

bamregion2GRanges <- function(bamfile, bamindex = bamfile, region = NULL, pairedEndReads = FALSE, min.mapq = 10, filterAltAlign = TRUE) {

  ## Check if bamindex exists
  bamindex.raw <- sub("\\.bai$", "", bamindex)
  bamindex <- paste0(bamindex.raw, ".bai")
  if (!file.exists(bamindex)) {
    bamindex.own <- Rsamtools::indexBam(bamfile)
    warning("Couldn't find BAM index-file ", bamindex, ". Creating our own file ", bamindex.own, " instead.")
    bamindex <- bamindex.own
  }

  ## Check if bam is truly paired ended in case pairedEndReads set to TRUE

  is.Paired <- Rsamtools::testPairedEndBam(file = bamfile, index = bamindex)

  if (pairedEndReads) {
    if (!is.Paired) {
      warning("You are trying to process single-ended BAM as paired-ended, Please set proper BAM directioanlity!!!")
    }
  } else {
    if (is.Paired) {
      warning("You are trying to process paired-ended BAM as single-ended, Please set proper BAM directioanlity!!!")
    }
  }

  ## read in reads data
  if (pairedEndReads) {
    print(Rsamtools::ScanBamParam(tag = c("XA", "HP"), which = range(region), what = c("mapq", "flag")))
    # suppressWarnings( data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(region), what=c('seq', 'qual','mapq','cigar'), flag=scanBamFlag(isDuplicate=F))) )
    data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index = bamindex, param = Rsamtools::ScanBamParam(tag = c("XA", "HP"), which = range(region), what = c("mapq", "flag")))
    print(data.raw)
  } else {
    data.raw <- GenomicAlignments::readGAlignments(bamfile, index = bamindex, param = Rsamtools::ScanBamParam(tag = c("XA", "HP"), which = range(region), what = c("mapq"), flag = scanBamFlag(isDuplicate = F)))
  }

  ## Second mate of the pair will inherit directionality from the first mate of the pair
  if (pairedEndReads) {
    data.first <- as(GenomicAlignments::first(data.raw), "GRanges")
    data.last <- as(GenomicAlignments::last(data.raw), "GRanges")
    strand(data.last) <- strand(data.first)
    data <- GenomicRanges::sort(c(data.first, data.last), ignore.strand = TRUE)
  } else {
    data <- as(data.raw, "GRanges")
  }

  ## Filter duplicates for pairedEndReads
  if (pairedEndReads) {
    bit.flag <- bitwAnd(1024, data$flag)
    mask <- bit.flag == 0
    data <- data[mask]
  }

  ## Filter by mapping quality
  if (!is.null(min.mapq)) {
    if (any(is.na(mcols(data)$mapq))) {
      warning(paste0(file, ": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
      mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
    }
    data <- data[mcols(data)$mapq >= min.mapq]
  }

  ## filter XA tag
  if (filterAltAlign) {
    data <- data[is.na(mcols(data)$XA)]
  }

  data <- GenomeInfoDb::keepSeqlevels(data, seqlevels(region), pruning.mode = "coarse")

  return(data)
}