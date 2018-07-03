#functions:
#' This function computes the dispersion parameters for W and C read counts in each segment
#'
#' @param prob.tab a table containing all segments with their states in all cells
#' , combined with all possible haplotypes and the corresponding C and W copy number
#' @param alpha The fraction of background reads in WW and CC states
#' @author Maryam Ghareghani
#' @export

add_dispPar <- function(prob.tab, alpha = 0.05)
{
  # set the dispersion parameters
  probs <- prob.tab
  probs[, disp_w := scalar * expected * nb_p / (1-nb_p) * Wcn * 0.5]
  probs[, disp_c := scalar * expected * nb_p / (1-nb_p) * Ccn * 0.5]
  
  # setting both W and C dispersion params for CN0 to alpha
  cn0_ridx <- which(probs$disp_c==0 & probs$disp_w==0)
  probs[cn0_ridx, disp_w := scalar * expected * nb_p / (1-nb_p) * alpha * 0.5]
  probs[cn0_ridx, disp_c := scalar * expected * nb_p / (1-nb_p) * alpha * 0.5]
  
  # rescaling W and C dispersion params for Wcn=0 cases with the parameter alpha
  Wcn0_ridx <- which(probs$disp_w==0)
  probs[Wcn0_ridx, disp_w:=disp_c*alpha]
  probs[Wcn0_ridx, disp_c:=disp_c*(1-alpha)]
  
  # rescaling W and C dispersion params for Ccn=0 cases with the parameter alpha
  Ccn0_ridx <- which(probs$disp_c==0)
  probs[Ccn0_ridx, disp_c:=disp_w*alpha]
  probs[Ccn0_ridx, disp_w:=disp_w*(1-alpha)]
  
  return(probs)
}


#' Compute the segment type given the segment strand state and the segment status
#' Returns W and C copy numbers, respectively
#' 
#' @param cellType The (majority) type of a cell that can have one these possible values: "ww","cc","wc","cw", or "?".
#' @param status A \code{vector} of length 4 containing {CN in hap1, inv CN in hap1, CN in hap2, inv CN in hap2} respectively.
#' @author Maryam Ghareghani
#' @export
#' 

getSegType = function(cellType, status)
{
  Wstatus = NULL
  if (cellType == "?")
  {
    segType = "?"
  }
  else
  {
    for (i in 1:2)
    {
      if (substr(cellType, i, i) == "W")
        Wstatus = paste0(Wstatus, "10")
      else
        Wstatus = paste0(Wstatus, "01")
    }
    Nw = Nc = 0
    
    for (i in 1:nchar(status))
    {
      Nw = Nw + as.integer(substr(status, i, i))*as.integer(substr(Wstatus, i, i))
      Nc = Nc + as.integer(substr(status, i, i))*(1-as.integer(substr(Wstatus, i, i)))
    }
    
    #segType = paste0(stringr::str_dup("W",Nw),stringr::str_dup("C",Nc))
    segType = c(Nw, Nc)
  }
  
  segType
}
