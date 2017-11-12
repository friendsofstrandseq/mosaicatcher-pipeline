#' Estimates p parameter and saves the empirical mean-var plot in the \code{directory}.
#' 
#' @inheritParams read.bams
#' @inheritParams nonzero.cov.bins 
#' @author Maryam Ghareghani
#' @export
#' 

estimateP <- function(counts, directory = NULL)
{
  numCells = (ncol(counts[[1]])-3)/2
  avg = NULL
  s = NULL
  for (i in 1:length(counts))
  {
    for (j in 1:numCells)
    {
      readcount = counts[[i]][,2*j+2]+counts[[i]][,2*j+3] # chr i, cell j
      avg = c(avg, mean(readcount))
      s = c(s, var(readcount))
    }
  }
  
  p = sum(avg*avg)/sum(s*avg)
  
  if (!is.null(directory))
  { 
    grDevices::pdf(paste0(directory, "mean-var-plot.pdf"))
    plot(avg, s, xlab = "mean", ylab = "var")    
    lines(avg, avg/p, col = "red")
    #TODO consider using ggplot so all plots use the grid graphics.
    #avg.df <- data.frame(avg=avg, pos=1:length(avg)) 
    #ggplot(avg.df) + geom_line(aes(x=pos, y=avg), color="red"))

    dev.off()
  }
  
  p
}


#' Estimates r parameters
#' 
#' @inheritParams nonzero.cov.bins 
#' @inheritParams getCNprob
#' @author Maryam Ghareghani
#' @export
#' 

estimateR <- function(counts, p)
{
  numCells = (ncol(counts[[1]])-3)/2
  avg = NULL
  s = NULL
  dispersion = matrix(, nrow = length(counts), ncol = numCells)
  for (i in 1:length(counts))
  {
    for (j in 1:numCells)
    {
      avg = mean(counts[[i]][,2*j+2]+counts[[i]][,2*j+3])
      dispersion[i,j] = avg*p/(1-p)
    }
  }
  
  dispersion
}


