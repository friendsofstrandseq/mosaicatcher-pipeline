#' Changes the format of all files to the format used in the SV calling wrapper function.
#'
#' @param RCfile The name of the file containing bin read counts
#' @param outputDir The directory containing the output files.
#' @param bamNamesFile Outputs the names of the bam files in the order they are used in the tool
#' @author Maryam Ghareghani
#' @export

# outputDir = "/local/home/mgharegh/research/HDhackathon/simulatedData/"
# RCfile = "/local/home/mgharegh/research/HDhackathon/simulatedData/counts.cov5.vaf0.5.small.p0.3.txt.gz"
# stateFile = "/local/home/mgharegh/research/HDhackathon/simulatedData/sces.cov5.vaf0.5.small.p0.3.txt" # "D2Rfb.strand_states.txt"
# bamNamesFile = "bamNames.txt"
# infoFile = "/local/home/mgharegh/research/HDhackathon/simulatedData/D2Rfb.50kb_fixed.info"
# breakpointsFile = "/local/home/mgharegh/research/HDhackathon/simulatedData/breakpoints.cov5.vaf0.5.small.p0.3.txt"

changeRCformat = function(RCfile, outputDir, bamNamesFile = "bamNames.txt")
{
  counts = data.table::fread(paste("zcat", RCfile))
  # newFormat is the table with all w counts first and then all the c counts
  newFormat = data.table::dcast(counts, chrom + start + end ~ cell, value.var = c("w", "c"))
  oldFormat = colnames(newFormat)[4:length(colnames(newFormat))]
  
  # exclude the extra chromosomes
  newFormat <- newFormat[grepl('^chr[0-9XY][0-9]?$', newFormat$chrom),]
  numcells = (ncol(newFormat)-3)/2
  # get the order of the columns to have W and C counts of the single cells together
  ord = NULL
  for (i in 1:numcells) (ord = c(ord, i, i+numcells))
  ord = c(1:3, ord + 3)
  data.table::setcolorder(newFormat, ord)
  
  # subset only autosomes
  newFormat = newFormat[which(sapply(newFormat$chrom, chrNumber) < 23),]
  
  # order the chromosomes numerically
  rowOrder = order(sapply(newFormat$chrom, chrNumber))
  newFormat = newFormat[rowOrder,]
  
  #naming
  colnames(newFormat)[1] = "chromosome"
  colnames(newFormat)[4:ncol(newFormat)] = paste0(rep(c("W", "C"), numcells), ceiling(1:(numcells*2)/2))
  
  # write mapping of cell name to cell number into file "bamNamesFile"
  oldColumOrder = data.frame(cell_id = 1:numcells,
                             cell_name = substr(oldFormat[1:numCells],3,nchar(oldFormat[1:numCells])))
  write.table(oldColumOrder, file = paste0(outputDir, bamNamesFile), quote = F, sep = "\t", row.names = F)
  
  newFormat
}

#' Changes the format of the cell types file and gives as output the cell types matrix
#'
#' @param stateFile The name of the file containing the cell types
#' @author Maryam Ghareghani
#' @export

changeCellTypesFormat = function(stateFile)
{
  #cellType = read.table(cellType, stringsAsFactors = F, header = T)
  d = data.table::fread(stateFile)
  d = unique(d)
  d = merge(d, d[, .(chrom_start = min(start), chrom_end = max(end)), by = chrom], by = "chrom")
  #d %>% group_by(chrom) %>% summarize(chrom_start = min(start), chrom_end = max(end))
  #merge(d, d %>% ..., by= "chrom")
  d = d[start == chrom_start & end == chrom_end,]
  x = data.table::dcast(d, chrom + start + end ~ sample + cell, value.var = "class")
  
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
#' @author Maryam Ghareghani
#' @export

changeNBparamsFormat = function(infoFile, K)
{
  info = data.table::fread(infoFile)
  p = info$nb_p[1]
  r = info$nb_r
  numcells = length(r)
  r = matrix(rep(r,K), ncol = numcells, byrow = T)
  
  list(p,r)
}


#' outputs the segment counts
#'
#' @param binRC bin read counts splitted by chromosomes
#' @param breakpointsFile The name of the breakpoint file
#' @param K The number of chromosomes (autosomes).
#' @param bin.size The size of the bins.
#' @author Maryam Ghareghani
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


# testing wether the names of the single cells have the same order in the RC file and the cell types file
# cellNames = as.numeric(sapply(colnames(newFormat)[4:(numcells+3)], substr, start = 3, stop = sapply(colnames(newFormat)[4:(numcells+3)], nchar)))
# cellNames2 = as.numeric(sapply(colnames(x)[4:(numcells+3)], substr, start = 7, stop = sapply(colnames(x)[4:(numcells+3)], nchar))) 
# table(cellNames == cellNames2)
# table(cellNames == info$cell)
# the order is the same
