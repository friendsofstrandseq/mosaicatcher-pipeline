#' Main function for the \pkg{\link{MaRyam}} package.
#'
#' Wrapper function to run SV calling after having the NB parmeters and segment read counts.
#'
#' @param bin.size the size of the bins
#' @param K The number of chromosomes (autosomes). #TODO I suggest that user can define what chromosome wants to analyze. So this parameter should take vector of chromosome IDs => c(1:22,'X') or paste0('chr', c(1:22,'X'))
#' @param maximumCN Maximum CN in the segments for SV calling. #TODO make sure parameters correspond to maxCN paramter in getCNprob
#' @param segmentsCounts The dataframe containing the segments chr and coordinates and W and C read counts over all cells
#' @param r The matrix of dispersion parameters of chromosomes(rows) and cells(columns) for CN=2 for a segment of size bin.size
#' @param p The p parameter of the NB distribution
#' @param cellTypes The matrix of cell types of chromosomes(rows) and cells(columns) TODO ... In function getSegType there is a parameter cellType. Are these paramters expecting the same thing? If not please use bit more informative name.
#' @param outputDir The outputDir containing the input and output files.
#' @param haplotypeMode a binary parameter showing wether we have the haplotype of WC cell types TODO ... #There is a haplotypeMode parameter in newSVcalling function. I guess they point to the same parameter.
#' @author Maryam Ghareghani
#' @export

#TODO If this is the main function please make sure that all relevant user-defined parameters can be set here!!!  

# non cmd
# outputDir = "/home/mgharegh/research/data/strand-seq/allCells-server/clean/HGSVC/HG00512/test/"
# segCountFile = "readCounts_5simple_inversions_HGSVC.txt"
# 
# # set parameters
# bin.size = 100000
# K = 22 # number of chromosomes
# maximumCN = 5

SVcalling_wrapperFunc = function(bin.size, K, maximumCN, segmentsCounts, r, p, cellTypes, outputDir, haplotypeMode)
{
  # blacklisting chr X and Y
  #segmentsCounts = segmentsCounts[which(sapply(segmentsCounts$chromosome, chrNumber) < 23),]
  
  colnames(cellTypes) = colnames(r)
  numCells = ncol(cellTypes)
  
  #the following line has bugs
  
  hapProbTables = list()
  normInvCNstatus = list()
  
  # store the subset of non SCE cells in each chromosome
  nonSCEcells = list()
  for (i in 1:K)
  {
    nonSCEcells[[i]] = which(cellTypes[i,] != "?")
  }
  
  hapStatus = NULL
  for (j in 0:maximumCN)
  {
    hapStatus = c(hapStatus, allStatus(3,j))
  }
  for (j in 1:length(hapStatus))
  {
    hapStatus[j] = paste(decodeStatus(hapStatus[j]), collapse = '')
  }
  
  aggProbTable = matrix(, nrow = nrow(segmentsCounts), ncol = length(hapStatus))
  colnames(aggProbTable) = hapStatus
  
  filterSeg <- NULL
  
  for (i in 1:nrow(segmentsCounts))
  {
    segCounts = segmentsCounts[i,]
    chr = chrNumber(segCounts[1,1])
    CN = getPossibleCNs(segCounts, p, as.numeric(r[chr,]), bin.size)
    print(paste("computing probabilities for", segCounts$chromosome, segCounts$start, segCounts$end))
    # filtering out all segments with weird CN (CN > maxCN)
    if (length(CN) > 0 && CN[1] <= maximumCN)
    {
      filterSeg <- c(filterSeg, i)
      # computing haplotype probabilities
      hapProbTables[[i]] = newgetCellStatProbabilities(hapStatus, counts = segCounts, chrCellTypes = as.character(cellTypes[chr,]), p, chrCellsDispPars = as.numeric(r[chr,]), binLength = bin.size, alpha = 0.05, haplotypeMode = haplotypeMode)
      
      #normalizing hapProbTable
      hapProbTables[[i]] = normalizeProbTable(hapProbTables[[i]])
      # regularization
      hapProbTables[[i]] = regularizeProbTable(hapProbTables[[i]])
      # computing aggregate Probabilities (based on non SCE cells only)
      for (h in 1:length(hapStatus))
      {
        aggProbTable[i,h] = sum(log(hapProbTables[[i]][h,nonSCEcells[[chr]]]))
      }
            
      #normalInvCNstatus
      normInvCNstatus[[i]] = getGenotypeProbTable(hapProbTables[[i]])
    }
  }
  
  # filter out small segments
  filterSeg = intersect(filterSeg, which(as.numeric(segmentsCounts$end) - as.numeric(segmentsCounts$start) > 100))
  

  # sort cells based on type in each chr
  chrOrder = list()
  for (i in 1:nrow(cellTypes))
  {
    chrOrder[[i]] = c(which(cellTypes[i,] == "cw"), which(cellTypes[i,] == "wc"), which(cellTypes[i,] == "ww"), which(cellTypes[i,] == "cc"))
  }
  
  # filter out all segments in the chromosomes that have SCE in all cells (If there is any) and send a message
  allSCEchrs <- match(TRUE, sapply(chrOrder, isEmpty))
  if (!is.na(allSCEchrs))
  {
    message(paste("Warning: chrs", allSCEchrs, "have SCE in all cells"))
    # keep only the segments that are in allSCEchrs
    newFilter <- which(is.na(match(sapply(segmentsCounts$chromosome, chrNumber), allSCEchrs)))
    filterSeg <- intersect(filterSeg, newFilter)
  }
  
  aggProbDF = data.frame()

  print("Making the dataframe for the heatmap...")
  for (i in filterSeg)
  {
    #print(i)
    if (is.null(hapProbTables[[i]])) # CN >= maxCN
    {
      next()
    }
    chr = chrNumber(segmentsCounts[i,1])
    df = data.frame(cell = chrOrder[[chr]], type = as.character(cellTypes[chr, chrOrder[[chr]]]), stringsAsFactors = F)
    
    Wcount = Ccount = NULL
    for (j in chrOrder[[chr]])
    {
      Wcount = c(Wcount, as.integer(segmentsCounts[i,2*j+2]))
      Ccount = c(Ccount, as.integer(segmentsCounts[i,2*j+3]))
    }
    
    df = cbind(df, Wcount)
    df = cbind(df, Ccount)
    
    cellsCNprob = matrix(, nrow = maximumCN+1, ncol = ncol(cellTypes))
    cellsCNprob[1,] = normInvCNstatus[[i]][1,]
    j = 2
    for (l in 1:maximumCN)
    {
      cellsCNprob[l+1,] = colSums(normInvCNstatus[[i]][j:(j+l),])
      j = j+l+1
    }
    
    for (j in 1:(maximumCN+1))
    {
      df = cbind(df, cellsCNprob[j,chrOrder[[chr]]])
    }
    
    for (h in 1:nrow(hapProbTables[[i]]))
    {
      df = cbind(df, hapProbTables[[i]][h,chrOrder[[chr]]])
    }
    
    aggProbVec = as.data.frame(c(segmentsCounts[i,1:3], ncol(cellTypes) + 1, "all", sum(df[,3]), sum(df[,4]), aggProbTable[i,]))
    names(aggProbVec) = c("chr", "start", "end", "cells", "types", "Wcount", "Ccount", hapStatus)
    aggProbDF = rbind(aggProbDF, aggProbVec)
    #print(paste("nrow =", nrow(aggProbDF)))
    
    segNames = segmentsCounts[i,1:3][rep(seq_len(nrow(segmentsCounts[i,1:3])), each=nrow(df)),]
    rownames(segNames) = NULL
    df = cbind(segNames, df)
    
    names(df) = c("chr", "start", "end", "cells", "types", "Wcount", "Ccount", paste0(rep("CN",maximumCN+1), as.character(0:maximumCN)), rownames(hapProbTables[[i]]))
    
    if (i == 1)
    {
      probTables = df
    } else
    {
      probTables = rbind(probTables, df)
    }
  }
  
  GTprobDF = getGenotypeProbDataFrame(probTables)

  print("writing the probabilities to the output files...")
  write.table(probTables, file = paste0(outputDir,"allSegCellProbs.table"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(GTprobDF, file = paste0(outputDir,"allSegCellGTprobs.table"), sep = "\t", quote = FALSE, row.names = FALSE)
  SVs = newSVcalling(aggProbDF)
  write.table(SVs, file = paste0(outputDir, "allSegSV.table"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(aggProbDF, file = paste0(outputDir,"allSegAggProbs.table"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# calling the function:
#SVcalling.wrapper.func(bin.size = 100000, K = 22, maximumCN = 5, outputDir = "/home/mgharegh/research/data/strand-seq/allCells-server/clean/HGSVC/HG00512/test/", segCountFile = "readCounts_5simple_inversions_HGSVC.txt")
#get segmentCounts as input instead of the segCountsFile
