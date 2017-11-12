#' Takes as an input an aggregate prob table and call each segment by the haplotype of the highest probability.
#' 
#' @param aggProbTable A \code{matrix} containing the logarithm of aggregate probabilities of haplotype states in segments, (rows correspond to segments and columns correspond to haplotype state).
#' @param haplotypeMode TODO ...
#' @author Maryam Ghareghani
#' @export

newSVcalling = function(aggProbTable, haplotypeMode = FALSE)
{
  SV = NULL
  names = colnames(aggProbTable)
  
  for (i in 1:nrow(aggProbTable))
  {
    mStat = which(as.numeric(aggProbTable[i,8:ncol(aggProbTable)]) == max(as.numeric(aggProbTable[i,8:ncol(aggProbTable)])))
    
    if (length(mStat) == 1)
    {
      SV = c(SV, names[7+mStat])
    } else if (!haplotypeMode & length(mStat) == 2 & names[7+mStat[1]] == sisterHaplotype(names[7+mStat[2]]))
    {
      SV = c(SV, names[7+mStat[1]])
    } else
    {
      SV = c(SV, NA)
    }
  }
  
  cbind(aggProbTable[,1:7], data.frame(invCNstat = SV))
}
