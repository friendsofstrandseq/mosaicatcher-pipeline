#' Read bam files from a directory and output them (the first mates only) as a list of \code{\link[GenomicRanges]{GRanges}} objects. It also writes the name of the bam files in a file.
#'
#' @param directory Directory containing the bam files.
#' @param bamFilenames Name of a file to write the names of the bam files in.
#' @param unq Only the unique alignments are given as output, if unq = TRUE.
#' @author Maryam Ghareghani
#' @export
#' 

readBams <- function(directory, bamFilenames = "", unq = TRUE)
{
  setwd(directory) #TODO !!! 
  files = list.files(pattern = "\\.bam$")  
  
  #It's not a good practice to change directory within the function. It's also not needed. See the code below. It will list bam files given the submitted path to BAM files. Please change. 		
  #TODO files <- list.files(directory, recursive=FALSE, pattern='.bam$') #consider to set full.names=TRUE

  if (bamFilenames != "")
  {
    write(noquote(files), file = bamFilenames, append = FALSE, sep = "\n")
  }
  
  grlist = list()
  for (i in 1:length(files))
  {
    print(paste("reading", files[i]))
    start.time = Sys.time()
    suppressWarnings(data.raw <- GenomicAlignments::readGAlignmentPairs(file = files[i],
                                                                    param=Rsamtools::ScanBamParam(what='mapq', flag=scanBamFlag(isDuplicate=F))))
                                                                   # , which = GRanges(paste0("chr",chrNum), IRanges(1,chrLen)))))
    print(paste("time for reading bam file:",Sys.time()-start.time))
    data.first <- as(GenomicAlignments::first(data.raw), 'GRanges')
    print(paste("time for subsetting the first mates:",Sys.time()-start.time))
    # quality filter--- getting unique reads
    if (unq)
    {
      data.first <- data.first[mcols(data.first)$mapq > 10]
    }
    print(paste("time for quality filtering:",Sys.time()-start.time))
    grlist[[i]] <- data.first
  }
  grlist
}
#TODO 1 need to be changed... We shouldn't count first pairs only
#TODO 2 Check the BAM import function below and implement in into your code. BAM import should work for single and paired end reads and should be able to take chromosomes as an input.

#' Import BAM file into GRanges
#'
#' Import aligned reads from a BAM file into a \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @param file Bamfile with aligned reads.
#' @param bamindex Bam-index file with or without the .bai ending. If this file does not exist it will be created and a warning is issued.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param min.mapq Minimum mapping quality when importing from BAM files.
#' @importFrom Rsamtools indexBam scanBamHeader ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments first last
#' @author David Porubsky
#' @export

bam2ranges <- function(file, bamindex=file, min.mapq=10, pairedEndReads=FALSE, chromosomes=NULL) {
  ## Check if bamindex exists
  bamindex.raw <- sub('\\.bai$', '', bamindex)
  bamindex <- paste0(bamindex.raw,'.bai')
  if (!file.exists(bamindex)) {
    bamindex.own <- Rsamtools::indexBam(file)
    warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
    bamindex <- bamindex.own
  }
  file.header <- Rsamtools::scanBamHeader(file)[[1]]
  chrom.lengths <- file.header$targets
  chroms.in.data <- names(chrom.lengths)
  if (is.null(chromosomes)) {
    chromosomes <- chroms.in.data
  }
  chroms2use <- intersect(chromosomes, chroms.in.data)
  if (length(chroms2use)==0) {
    chrstring <- paste0(chromosomes, collapse=', ')
    stop('The specified chromosomes ', chrstring, ' do not exist in the data. Please try ', paste(paste0('chr',chromosomes), collapse=', '), ' instead.')
  }
  ## Issue warning for non-existent chromosomes
  diff <- setdiff(chromosomes, chroms.in.data)
  if (length(diff)>0) {
    diffs <- paste0(diff, collapse=', ')
    warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data.'))
  }
  ## Import the file into GRanges
  gr <- GenomicRanges::GRanges(seqnames=chroms2use, ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))

  if (pairedEndReads) {
      data.raw <- GenomicAlignments::readGAlignmentPairs(file, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))
  } else {
      data.raw <- GenomicAlignments::readGAlignments(file, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))
  }
   
  ## Second mate of the pair will inherit directionality from the first mate of the pair
  if (pairedEndReads) {
    data.first <- as(GenomicAlignments::first(data.raw), 'GRanges')
    data.last <- as(GenomicAlignments::last(data.raw), 'GRanges')
    strand(data.last) <- strand(data.first)
    data <- sort(c(data.first, data.last))
  } else {
    data <- as(data.raw, 'GRanges')
  }
  
  ## Filter by mapping quality
  if (!is.null(min.mapq)) {
    if (any(is.na(mcols(data)$mapq))) {
      warning(paste0(file,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
      mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
    }
    data <- data[mcols(data)$mapq >= min.mapq]
  }
  seqlevels(data) <- seqlevels(gr)
  return(data)
}
