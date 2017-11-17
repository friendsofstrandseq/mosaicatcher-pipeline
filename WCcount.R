#' Count the number of Watson and Crick reads and return a \code{dataframe} of read counts.
#' 
#' @param cells.alignmemts A \code{list} of \code{\link[GenomicRanges]{GRanges}} objects with Strand-specific read data. #TODO Is the object a GRangesList or a normal list of GRanges??? If is is a GRangesList please use \code{\link[GenomicRanges]{GRangesList}}
#' @param segments A list \code{\link[GenomicRanges]{GRanges}} objects (genomic intervals).
#' @author Maryam Ghareghani
#' @export
#' 

WCreadCounts <- function(segments, cells.alignments)
{
  numCells = length(cells.alignments)
  
  df = data.frame(chromosome = seqnames(segments), start = start(segments), end = end(segments))
  
  for (i in 1:numCells)
  {
    reads <- cells.alignments[[i]]
    wdata <- reads[strand(reads) == "-"]
    cdata <- reads[strand(reads) == "+"]
    
    df[[paste0("W",i)]] = GenomicRanges::countOverlaps(segments, wdata)
    df[[paste0("C",i)]] = GenomicRanges::countOverlaps(segments, cdata)
  }
  df
}


#' Splits a \code{data.frame} of read counts based on chromosomes and returns a \code{list} containing the splited dataframes.
#' 
#' @inheritParams nonzeroCovBins
#' @author Maryam Ghareghani
#' @export
#'

#TODO Not sure how the counts object looks like but I think this can be done using: split(counts, <chromosome>); Split returns list of data.frames splited based on chromosome ID. Please check.  

splitChromosomes <- function(counts)
{
  chr = unique(counts$chromosome)
  list.counts = list()
  for (i in 1:length(chr))
  {
    list.counts[[i]] = counts[counts$chromosome == as.character(chr)[i],]
  }
  list.counts
}

# considers only autosomes
