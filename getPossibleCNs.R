#' Returns all CNs with maximum probability.
#' 
#' @inheritParams getCNprob
#' @param segCounts W and C counts of one segment (a count \code{data.frame} with one row).
#' @param chrCellsDispPars The dispersion parameters of cells for the chromosome where the segment resides.
#' @author Maryam Ghareghani
#' @export
#' 

getPossibleCNs = function(segCounts, p, chrCellsDispPars, binLen, alpha)
{
  CN = NULL
  totalR = sum(chrCellsDispPars)
  totalCount = sum(as.integer(segCounts[,4:ncol(segCounts)]))
  segLen = as.integer(segCounts[,3]) - as.integer(segCounts[,2]) + 1
  
  CNp = getCNprob(totalCount, p, totalR, segLen, binLen)
  
  if (sum(CNp) > 0) # CN <= maxCN = 5
  {
    CN = which(CNp == max(CNp))-1
  }
  
  CN
}


#' Generate a CN probability vector.
#' 
#' @param readCount Total number of read counts (W+C) in all cells in this segment.
#' @param p The p parameter of the NB distribution.
#' @param r Dispersion parameter for the total number of read counts (W+C in all cells) in a bin with a diploid copy number (2).
#' @param segLen The length of the segment.
#' @param binLen The length of the bins used for estimating the NB parameters.
#' @param maxCN The maximum number of CNs.
#' @param alpha The coefficient of NB dispersion parameter for the number of background reads.
#' @author Maryam Ghareghani
#' @export
#' 

#TODO it seems that maxCN is set as default here and user cannot change it via wrapper function. Also function SVcalling.wrapper.func takes a parameter maximumCN which seems to be related. Please synchronize this.

getCNprob = function(readCount, p, r, segLen, binLen, maxCN = 5, alpha = 0.05)
{
  CNprob = NULL
  disp = (r/2)*(segLen/binLen) # disp = dispersion par for copy number 1 for this segLen
  
  for (i in 0:maxCN)
  {
    CNprob = c(CNprob, dnbinom(readCount, size = disp*(max(i,alpha)), prob = p))
  }
  
  if (sum(CNprob) != 0)
  {
    CNprob = CNprob/sum(CNprob) # normalization
  }
  else
  {
    expRC = (1-p)*disp/p
    CN = round(readCount/expRC)
    CNprob = c(rep(0,CN),1)
  }
  
  CNprob
}
