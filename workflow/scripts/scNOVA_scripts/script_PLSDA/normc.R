# library(pracma)
#' Normaliz the columns of x to a length of 1.
#'
#' @param x n*p matrix
#'
#' @return xn normalized result
#' @export normc
#'
#' @examples
#' #ex1.
#' m <- matrix(1:4,2,2,byrow=TRUE)
#' normc(m)
#' #ex2.
#' n <- matrix(rnorm(100,10,1),10,10)
#' normc(n)
normc <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  xn <- x
  for(j in 1:p){
    tempscale <- norm(matrix(x[,j]),type="f")
    xn[,j] <- matrix(x[,j])/tempscale
  }
  return(xn)
}
