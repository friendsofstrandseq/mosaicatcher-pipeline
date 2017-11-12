#' Filters the list of data.frames based on the list od indexes.
#' 
#' @param list.df A \code{list} of data.frames.
#' @param list.idx A \code{list} of indexes for subsetting rows of the \code{data.frame}.
#' @author Maryam Ghareghani
#' @export
#' 

filt <- function(list.df, list.idx)
{
  new.list.df = list()
  for (i in 1:length(list.df))
  {
    new.list.df[[i]] = list.df[[i]][list.idx[[i]],]
  }
  new.list.df
}
#TODO need to be changed... We shouldn't count first pairs only


#' Outputs the list of indices of the bins with non zero coverage.
#' 
#' @param counts A \code{data.frame} containing W and C read counts in some segments.
#' @author Maryam Ghareghani
#' @export
#' 

nonzero.cov.bins <- function(counts)
{
  idx = list()
  for (i in 1:length(counts))
  {
    idx[[i]] = which(rowSums(counts[[i]][,4:ncol(counts[[i]])]) != 0)
  }
  idx
}
#TODO need to be changed... We shouldn't count first pairs only


#' Outputs the list of indices of the bins having the fraction of low mapq reads less than alpha.
#'
#' @inheritParams nonzero.cov.bins
#' @inheritParams getCNprob
#' @param uniqe.counts A \code{data.frame} containing W and C unique read counts in segments.
#' @author Maryam Ghareghani
#' @export
#' 

unique.mappable.bins <- function(counts, unique.counts, alpha = 0.05)
{
  idx = list()
  for (i in 1:length(counts))
  {
    idx[[i]] = which(rowSums(unique.counts[[i]][,4:ncol(unique.counts[[i]])])/rowSums(counts[[i]][,4:ncol(counts[[i]])]) > (1-alpha))
  }
  idx
}
#TODO need to be changed... We shouldn't count first pairs only
