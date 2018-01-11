#' Changes the format of all files to the format used in the SV calling wrapper function.
#'
#' @param RCfile The name of the file containing bin read counts
#' @param outputDir The directory containing the output files.
#' @param bamNamesFile Outputs the names of the bam files in the order they are used in the tool
#' @author Maryam Ghareghani, Sascha Meiers
#' @export

changeRCformat = function(RCfile, outputDir, bamNamesFile = "bamNames.txt")
{
  counts = data.table::fread(paste("zcat", RCfile))
  # newFormat is the table with all w counts first and then all the c counts
  newFormat = data.table::dcast(counts, chrom + start + end ~ cell, value.var = c("w", "c"))
  oldFormat = colnames(newFormat)[4:length(colnames(newFormat))]
  
  # exclude the extra chromosomes
  newFormat <- newFormat[grepl('^chr[0-9XY][0-9]?$', newFormat$chrom),]
  numCells = (ncol(newFormat)-3)/2
  # get the order of the columns to have W and C counts of the single cells together
  ord = NULL
  for (i in 1:numCells) (ord = c(ord, i, i+numCells))
  ord = c(1:3, ord + 3)
  data.table::setcolorder(newFormat, ord)
  
  # subset only autosomes
  newFormat = newFormat[which(sapply(newFormat$chrom, chrNumber) < 23),]
  
  # order the chromosomes numerically
  rowOrder = order(sapply(newFormat$chrom, chrNumber))
  newFormat = newFormat[rowOrder,]
  
  #naming
  colnames(newFormat)[1] = "chromosome"
  colnames(newFormat)[4:ncol(newFormat)] = paste0(rep(c("W", "C"), numCells), ceiling(1:(numCells*2)/2))
  
  # write mapping of cell name to cell number into file "bamNamesFile"
  oldColumOrder = data.frame(cell_id = 1:numCells,
                             cell_name = substr(oldFormat[1:numCells],3,nchar(oldFormat[1:numCells])))
  write.table(oldColumOrder, file = paste0(outputDir, bamNamesFile), quote = F, sep = "\t", row.names = F)
  
  list(binRC=newFormat, cellNames = as.character(oldColumOrder$cell_name))
}

#' Changes the format of the cell types file and gives as output the cell types matrix
#'
#' @param stateFile The name of the file containing the cell types
#' @param cellNames The names of the single cells in the correct order
#' @author Maryam Ghareghani, Sascha Meiers
#' @export

changeCellTypesFormat = function(stateFile, cellNames)
{
  #cellType = read.table(cellType, stringsAsFactors = F, header = T)
  d = data.table::fread(stateFile)
  d = unique(d)
  # adding two columns to the end including the start and end of chromosomes
  d = merge(d, d[, .(chrom_start = min(start), chrom_end = max(end)), by = chrom], by = "chrom")
  # kick out the SCE cells
  d = d[start == chrom_start & end == chrom_end,]
  x = data.table::dcast(d, chrom + start + end ~ cell, value.var = "class")
  
  # make sure if all cells are reported
  stopifnot(unique(sort(d$cell)) == unique(sort(cellNames)))
  
  # order the columns based on the order of the cellNames
  names <- colnames(x)[4:ncol(x)]
  m <- match(cellNames, names)
  x <- x[,c(1:3, m+3),with = F]
  
  
  # exclude the extra chromosomes
  x <- x[grepl('^chr[0-9XY][0-9]?$', x$chrom),]
  
  # sorting rows and exclude X and Y
  # order the chromosomes numerically
  rowOrder = order(sapply(x$chrom, chrNumber))
  # subset only autosomes
  x = x[rowOrder[rowOrder < 23],]
  
  cellTypes = as.matrix(x[,4:ncol(x)])
  cellTypes[which(is.na(cellTypes))] = "?"
  
  tolower(cellTypes)
}

#' Changes the format of the NB parameters and gives as output p parameter and a matrix containing r parameters
#'
#' @param infoFile The name of the file containing the NB parameters
#' @param K The number of chromosomes (autosomes).
#' @param cellNames The names of the single cells in the correct order
#' @author Maryam Ghareghani
#' @export

changeNBparamsFormat = function(infoFile, K, cellNames)
{
  info = data.table::fread(infoFile)
  p = info$nb_p[1]
  r = info$nb_r
  numCells = length(r)
  r = matrix(rep(r,K), ncol = numCells, byrow = T)
  
  list(p,r)
}


#' outputs the segment counts
#'
#' @param binRC bin read counts splitted by chromosomes
#' @param breakpointsFile The name of the breakpoint file
#' @param K The number of chromosomes (autosomes).
#' @param bin.size The size of the bins.
#' @author Maryam Ghareghani, Sascha Meiers
#' @export

getSegReadCounts = function(binRC, breakpointsFile, K, bin.size)
{
  seg = utils::read.table(breakpointsFile, stringsAsFactors = F, colClasses = c("integer", "character", "integer"), header = T)[,2:3]
  colnames(seg) = c("chromosome", "breakpoint")
  segRC = data.frame()
  
  for (k in 1:K)
  {
    print(paste("k =", k))
    chrSegs = seg[seg$chromosome == paste0("chr",k),]
    if (nrow(chrSegs) < 2)
      next()
    
    
    for (i in 2:nrow(chrSegs)) # assumption: there are at least two breakpoints in each chromosome
    {
      start.bin.idx = chrSegs$breakpoint[i-1] + 1
      end.bin.idx = chrSegs$breakpoint[i]
      
      df = data.frame(chromosome = paste0("chr",k), start = bin.size*(start.bin.idx), end = bin.size*(end.bin.idx+1),stringsAsFactors = F)
      subRows = (start.bin.idx:end.bin.idx)

      df = cbind(df, t(as.data.frame(colSums(as.matrix(binRC[[k]][subRows, 4:ncol(binRC[[k]])])))))
      rownames(df)  = NULL
      
      segRC = rbind(segRC, df)
    }
  }
  
  segRC
}
