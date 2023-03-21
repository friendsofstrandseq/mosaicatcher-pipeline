auto_R <- function(x){
  library(matrixStats)
  m <- nrow(x); n <- ncol(x);
  mx    = colMeans(x);
  stdx  = colSds(x);
  ax    = (x-t(matrix(1, ncol(x), nrow(x))*mx)) / t(matrix(1, ncol(x), nrow(x))*stdx);
  return(ax)
}