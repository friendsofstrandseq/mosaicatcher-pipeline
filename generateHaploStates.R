#' Returns the first binary string with n 01 and m 1s.
#' 
#' @param n The number of 0s.
#' @param m The number of 1s.
#' @author Maryam Ghareghani
#' @export
#' 

#TODO This function is very short. I would define this function at the beginning of functions defined below. 

initialState = function(n, m)
{
  paste0(stringr::str_dup("0",n), stringr::str_dup("1",m))
}


#' Computes the next binary status of a given status with equal number of 0s and 1s.
#' 
#' @param currentState The current binary status.
#' @author Maryam Ghareghani
#' @export
#' 

getNextState = function(currentState)#, n, m)
{
  nextState = currentState
  pos = stringr::str_locate_all(currentState, "01")[[1]][,1]
  if (length(pos) == 0)
  {
    return(FALSE) # can't be incremented
  }
  else
  {
    pos = pos[length(pos)] # last occurence of the pattern
    nextState = paste0(substr(nextState,1,pos-1),"10",substr(nextState,pos+2,nchar(nextState)))
    nextState = paste(nextState, collapse = "")
    
    if (pos + 1 < nchar(nextState))
    {
      c1 = stringr::str_count(substr(nextState, pos+2, nchar(nextState)),"1")
      c0 = nchar(nextState) - pos - 1 - c1
      nextState = paste0(substr(nextState,1,pos+1), initialState(c0,c1))
    }
  }
  nextState
}


#' Generates all binary strings with n 0s and m 1s.
#' 
#' @inheritParams initialState
#' @author Maryam Ghareghani
#' @export
#' 


allStatus = function(n,m)
{
  allStat = NULL
  status = initialState(n,m)
  
  while(status != FALSE)
  {
    allStat = c(allStat, status)
    status = getNextState(status)
  }
  
  allStat
}


#' Converts a binary status to an integer vector status.
#' 
#' @param binaryStatus A binary \code{vector} in which the number of 1s between every two consecutive zeros (or before the first zero or after the last zero) indicates a copy number.
#' @author Maryam Ghareghani
#' @export
#' 

decodeStatus = function(binaryStatus)
{
  status = NULL
  pos = stringr::str_locate_all(binaryStatus, "0")[[1]][,1]
  pos = c(0, pos, nchar(binaryStatus)+1)
  
  for (i in 1:length(pos)-1)
  {
    status = c(status, pos[i+1]-pos[i]-1)
  }
  status
}


#' Takes as an input a haplotype state and returns the opposite haplotype state which has the same genotype.
#' 
#' @param hapState A decoded haplotype state.
#' @author Maryam Ghareghani
#' @export
#' 

sisterHaplotype = function(hapState)
{
  sisterHap = ""
  for (i in 1:(nchar(hapState)/4))
  {
    sisterHap = paste0(sisterHap, substr(hapState,i+2,i+3), substr(hapState,i,i+1))
  }
  sisterHap
}


#TODO Remove or at least rename it later!!!
#' Computes the total normal CN and inverted CN status corresponding to the input decoded (non-binary) haplotype status.
#' 
#' Note that it doesn't compute the genotype (distribution of the CNs in the haplotypes).
#'
#' @param haplotyeStatus A decoded (non-binary) haplotype status.
#' @author Maryam Ghareghani
#' @export
#' 

getGenotypeStatus = function(haplotypeStatus)
{
  hapStatus = as.integer(strsplit(haplotypeStatus,"")[[1]])
  genotypeStatus = NULL
  for (i in 1:(length(hapStatus)/4))
  {
    genotypeStatus = c(genotypeStatus, hapStatus[4*i-3]+hapStatus[4*i-1])
    genotypeStatus = c(genotypeStatus, hapStatus[4*i-2]+hapStatus[4*i])
  }
  
  paste(genotypeStatus, collapse = "")
}
