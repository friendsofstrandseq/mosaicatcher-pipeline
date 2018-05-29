#' Estimate mappability in user defined bins.
#'
#' This function takes user defined bin sizes and estimates mappable fraction based on high coverage sequencing data.
#' 
#' @param bamFile A list of files that contains \code{\link{BreakPoint}} objects.
#' @param export A location to write the results (if 'exportBed' set to TRUE)
#' @param min.mapq A minimal mapping quality of a read to be considered.
#' @param pairedEndReads Set to TRUE if reads are paired.
#' @param bin.sizes A single value or a vector of bin sizes (in bp) to count reads in.
#' @param chunkSize A size of genomic region to load in order to prevent exceeding memory for large BAM files.
#' @param chromosomes A list of chromosomes to process.
#' @param exportBed Set to TRUE if you want to store binned data in a bed file (Set also 'export' location)
#' @return A \code{list} object of binned data for all submitted 'bin.sizes'.
#' @author David Porubsky

## Run
#estimateMappability(bamFile="/media/porubsky/5D7CBC85372E82B0/StrandSeqNation/RPE_WT_srt_dedup.bam", export="/media/porubsky/5D7CBC85372E82B0/StrandSeqNation/", min.mapq=10, pairedEndReads=FALSE, bin.sizes=c(50000,100000,200000), chunkSize=10000000, chromosomes = paste0('chr', c(1:22,"X")), exportBed=TRUE)

estimateMappability <- function(bamFile, export=".", min.mapq=10, pairedEndReads=FALSE, bin.sizes=50000, chunkSize=10000000, chromosomes=NULL, exportBed=TRUE) {
  
  ## Get chromosome length
  file.header <- Rsamtools::scanBamHeader(bamFile)[[1]]
  chrom.lengths <- file.header$targets
  chroms.in.data <- names(chrom.lengths)
  if (is.null(chromosomes)) {
    chromosomes <- chroms.in.data
  }
  chroms2use <- intersect(chromosomes, chroms.in.data)
  
  #process data in chunks to save memory
  chunks <- unlist(tileGenome(chrom.lengths[chroms2use], tilewidth = chunkSize))
  
  cumCov.perLib <- GRanges()
  for (i in 1:length(chunks)) {
    chunk <- chunks[i]
    
    data.raw <- GenomicAlignments::readGAlignments(bamFile, param=Rsamtools::ScanBamParam(which=range(chunk), what='mapq', flag=scanBamFlag(isDuplicate=FALSE)))
    data <- as(data.raw, 'GRanges')
    ## Filter by mapping quality
    if (!is.null(min.mapq)) {
      if (any(is.na(mcols(data)$mapq))) {
        warning(paste0(file,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
        mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
      }
      data <- data[mcols(data)$mapq >= min.mapq]
    }
    
    if (length(cumCov.perLib)==0) {
      cumCov.perLib <- reduce(data)
    } else {
      cumCov.perLib <- reduce(c(cumCov.perLib, data[,0]))
    }
  }
  cumCov.perLib <- keepSeqlevels(cumCov.perLib, chromosomes)
  cov <- coverage(cumCov.perLib)
  
  binned.binsizes <- list() 
  for (bin.len in bin.sizes) {
    message("Working on binsize ", bin.len, "bp")
    
    chrom.lengths.floor <- floor(chrom.lengths[chromosomes] / bin.len) * bin.len
    binned <- unlist(tileGenome(chrom.lengths.floor, tilewidth = bin.len), use.names = FALSE)
    
    binned.mappable <- GRangesList()
    seqlevels(binned.mappable) <- chromosomes
    for(chr in names(cov)) {
      chr.cov <- cov[[chr]]
      binned.chr <- binned[seqnames(binned) == chr]
      bin.cov <- Views(chr.cov, ranges(binned.chr))
      cov.bases <- viewApply(bin.cov, function(x) sum(as.vector(x > 0)))
      binned.chr$mappable.frac <- cov.bases / width(binned.chr)
      suppressWarnings( binned.mappable[[chr]] <- binned.chr )
    }
    binned.mappable <- unlist(binned.mappable, use.names = FALSE)
    binned.binsizes[[paste0(format(bin.len, scientific = FALSE), "bin")]] <- binned.mappable
    
    if (exportBed) {
      binned.mappable.df <- as(binned.mappable, 'data.frame')
      binned.mappable.bed <- binned.mappable.df[,c('seqnames', 'start', 'end', 'mappable.frac')]
      
      filename <- paste0(format(bin.len, scientific=FALSE),"bp_mappableFrac_WGS.bed")
      destination <- file.path(export, filename)
      write.table(binned.mappable.bed, file = destination, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    }
  }
  return(binned.binsizes)
}  


#####################
## HELPER FUNCTION ##
#####################

reformat <- function(x) {
  out_list <- list() 
  for ( i in seq(1, length(x), 2) ) {
    out_list[[i]] <- c(x[i], x[i+1])
  }
  mt <- do.call("rbind",out_list)
  df <- data.frame(mt)
  colnames(df) <- c("start", "end")
  df
}  
