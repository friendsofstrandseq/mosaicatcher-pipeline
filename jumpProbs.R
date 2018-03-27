#' Computes the probabilities of having the same genotype/haplotype for all single cells (jump probabilities)
#' returns a list containing the lists of probability tables of segments for each chromosome
#' In the output probability tables, the input and output jump probabilities are included for internal (not first and last) segments
#' 
#' @param probTable.l.chroms 
#' @author Maryam Ghareghani
#' @export
#' 

addJumpProbs <- function(probTable, simplified=TRUE)
{
  maximumCN <- length(which(startsWith(colnames(probTable), "CN")))-1
  n <- ncol(probTable)
  
  if (simplified)
  {
    s = 8
  } else {
    s = maximumCN+9
  }
  
  # split probTable by chrom
  probTable.l.chroms <- split(probTable, factor(probTable$chr, levels = unique(probTable$chr)))
  
  # split probTable.l.chrom by segment
  for (k in 1:length(probTable.l.chroms))
  {
    ID <- paste0(probTable.l.chroms[[k]]$start, "_", probTable.l.chroms[[k]]$end)
    ID <- factor(ID, levels=unique(ID))
    probTable.l.chroms[[k]] <- split(probTable.l.chroms[[k]], ID)
  }
  
  jump.probs <- list()
  for (k in 1:length(probTable.l.chroms)) {
    jump.probs[[k]] <- list()
    if (length(probTable.l.chroms[[k]]) == 1)
      next()
    for (i in 1:(length(probTable.l.chroms[[k]])-1)) {
      # computing the jump prob from segment i to i+1
      prod.probs <- probTable.l.chroms[[k]][[i]][,s:n] * probTable.l.chroms[[k]][[i+1]][,s:n]
      jump.probs[[k]][[i]] <- 1-rowSums(prod.probs)
      
      # adding jump probs to the prob table
      if (i > 1)
      {
        probTable.l.chroms[[k]][[i]] <- cbind(probTable.l.chroms[[k]][[i]], 
                data.frame(input_jump_p = jump.probs[[k]][[i-1]], output_jump_p = jump.probs[[k]][[i]]))
      }
    }
  }
  
  # return the list of segments prob table including the jump probabilities
  return(probTable.l.chroms)
}
