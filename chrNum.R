#' returns the number of a chromosome, given the input string
#' @param chr a string like "chri"
#' @author Maryam Ghareghani
#' @export
#' 

chrNumber = function(chr)
{
  c = 0
  chr = substr(chr, 4, nchar(chr)) # remove "chr" from the chr name
  if (chr == "X")
    c = 23
  else if (chr == "Y")
    c = 24
  else
    c = as.numeric(chr)
  c
}
