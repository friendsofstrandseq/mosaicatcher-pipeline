#' Counts the number of Watson and Crick reads of single cells in bins and segments and fitting NB distribution.
#'
#' @param segmentsFile The name of the segments bed file.
#' @param temDir The directory containing the chromosomes length file.
#' @param bamDir The directory containing all of the bam files.
#' @param directory The directory containing the input and output files.
#' @param bin.size The size of the bins.
#' @param K the number of chromosomes (autosomes).
#' @author Maryam Ghareghani
#' @export

#args = commandArgs(trailingOnly=TRUE)
#tempDir = args[1]
#bamDir = args[2]
#directory = args[3]

# set parameters
#bin.size = 100000
#K = 22 # number of chromosomes

countAndNBfit.wrapper.func = function(segmentsFile, tempDir, bamDir, directory, bin.size, K = 22)
{
  # get length of chromosomes
  chrLens = utils::read.table(paste0(tempDir,"chrLens.data"))[,1]
  
  # read bam files
  start.time = Sys.time()
  cell.alignments = read.bams(bamDir, paste0(directory,"bamFilenames.data"), FALSE)
  print(Sys.time() - start.time)
  numCells = length(cell.alignments)
  
  cell.unique.alignments = list()
  for (i in 1:length(cell.alignments))
  {
    cell.unique.alignments[[i]] = cell.alignments[[i]][mcols(cell.alignments[[i]])$mapq > 10]
  }
  print(Sys.time() - start.time)
  
  # Getting cell types
  cellTypes = utils::read.table(file = paste0(directory, "cellTypes.data"), stringsAsFactors = F)
  
  # make consecutive intervals of length bin.size
  bins = data.frame()
  for (i in 1:K)
  {
    numbins = floor(chrLens[i]/bin.size)
    df = as.data.frame(successiveIRanges(rep(bin.size, numbins)))
    df = cbind(chromosome = rep(paste0("chr",i), numbins), df)
    bins = rbind(bins, df)
  }
  
  # constructing a granges object from bins
  segments <- GenomicRanges::GRanges(seqnames=bins$chromosome, ranges=IRanges(start=bins$start, end=bins$end))
  
  # count the number of W and C reads in the bins
  counts = WCreadCounts(segments, cell.alignments)
  unique.counts = WCreadCounts(segments, cell.unique.alignments)
  
  # output bin read counts to a file
  utils::write.table(counts, file = paste0(directory,"binReadCounts.data"), quote = FALSE, row.names = FALSE)
  utils::write.table(unique.counts, file = paste0(directory,"binUniqeReadCounts.data"), quote = FALSE, row.names = FALSE)
  
  counts = split.chromosomes(counts)
  unique.counts = split.chromosomes(unique.counts)
  
  # blacklisting
  nonzeroIndex = nonzero.cov.bins(counts)
  counts = filt(counts, nonzeroIndex)
  unique.counts = filt(unique.counts, nonzeroIndex)
  
  mappable.bins.index = unique.mappable.bins(counts, unique.counts)
  counts = filt(counts, mappable.bins.index)
  unique.counts = filt(unique.counts, mappable.bins.index)
  
  # estimate NB parameters
  p = estimateP(unique.counts, directory)
  disp = estimateR(unique.counts, p)
  rownames(disp) = paste0(rep("chr",K), 1:K)
  colnames(disp) = paste0(rep("cell", numCells), 1:numCells)
  
  start.time = Sys.time()
  if (!dir.exists(paste0(directory,"NBfitPlots")))
  {
    system(paste0("mkdir ", directory,"NBfitPlots"))
  }
  NBfitplots(paste0(directory, "NBfitPlots/"), unique.counts, cellTypes, p, disp, bin.size)
  print(Sys.time()-start.time)
  
  # outputing the NB parameters
  write(p, file = paste0(directory, "p.data"))
  write.matrix(disp, paste0(directory, "r.data"))
  
  # counting W and C reads in segments
  df.segm = read.table(paste0(directory, segmentsFile), header=F)
  colnames(df.segm) = c("chromosome", "start", "end")
  segments <- GenomicRanges::GRanges(seqnames=df.segm$chromosome, ranges=IRanges(start=df.segm$start, end=df.segm$end))
  seg.unique.counts = WCreadCounts(segments, cell.unique.alignments)
  utils::write.table(seg.unique.counts, file = paste0(directory,"readCounts_",segmentsFile), quote = FALSE, row.names = FALSE)
}

#countAndNBfit.wrapper.func("segments.bed", tempDir = "/local/home/mgharegh/research/data/strand-seq/allCells-server/clean/", bamDir = "/local/home/mgharegh/research/data/strand-seq/allCells-server/clean/HGSVC/first2HG00512bams/", directory = "/local/home/mgharegh/research/data/strand-seq/allCells-server/clean/HGSVC/HG00512test/", bin.size = 100000)
