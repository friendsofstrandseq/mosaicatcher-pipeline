# parameters:
maximumCN <- 4
haplotypeMode <- F
regularizationParameter <- 1e-10

###-------------------------------------

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
  new.probs <- prob.tab
  new.probs[, disp_w := scalar*nb_r*Wcn*(to-from+1)*(0.5)]
  new.probs[, disp_c := scalar*nb_r*Ccn*(to-from+1)*(0.5)]
  
  # setting both W and C dispersion params for CN0 to alpha
  cn0_ridx <- which(new.probs$disp_c==0 & new.probs$disp_w==0)
  new.probs[cn0_ridx, disp_w:=scalar*nb_r*alpha*(to-from+1)]
  new.probs[cn0_ridx, disp_c:=scalar*nb_r*alpha*(to-from+1)]
  
  # rescaling W and C dispersion params for Wcn=0 cases with the parameter alpha
  Wcn0_ridx <- which(new.probs$disp_w==0)
  new.probs[Wcn0_ridx, disp_w:=disp_c*alpha*2]
  new.probs[Wcn0_ridx, disp_c:=disp_c*(1-alpha)*2]
  
  # rescaling W and C dispersion params for Ccn=0 cases with the parameter alpha
  Ccn0_ridx <- which(new.probs$disp_c==0)
  new.probs[Ccn0_ridx, disp_c:=disp_w*alpha*2]
  new.probs[Ccn0_ridx, disp_w:=disp_w*(1-alpha)*2]
  
  return(new.probs)
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
    
    segType = paste0(stringr::str_dup("W",Nw),stringr::str_dup("C",Nc))
  }
  
  #list(segType, Wcn = Nw, Ccn = Nc)
  c(Nw, Nc)
}
##-------------------------

# defining the vector of all possible haplotypes
hapStatus <- NULL
for (j in 0:maximumCN)
{
  hapStatus <- c(hapStatus, allStatus(3,j))
}
for (j in 1:length(hapStatus))
{
  hapStatus[j] <- paste(decodeStatus(hapStatus[j]), collapse = '')
}

# creating a datatable containing all possible combinations of strand states and haplotypes, and setting their segTypes
hapStrandStates <- data.table()
for (st in c("CC","WW","WC","CW"))
{
  hapStrandStates <- rbindlist(list(hapStrandStates, 
                     data.table(state=st, haplotype=hapStatus, 
                     segtype=t(sapply(hapStatus, function(x) getSegType(st, x))))))
}
colnames(hapStrandStates)[3:4]=c("Wcn", "Ccn")
setkey(hapStrandStates, state)

# I don't know how to merge these two datatables (probs and hapStrandStates)
#take the only solution that I have in mind

# kick out the segs with sces
probs <- probs[state!="sce"]
#probs <- probs[order(sample, cell, chrom, from, to)]
# compute the numer of segments with different strand states in probs
strand.count <- probs[,.N,by=state]
setkey(strand.count,state)
strand.count <- strand.count$N
# expand the hapStrandStates datatable
rep.rows <- NULL
for (i in 1:4)
{
  rep.rows <- c(rep.rows,
              rep((1:length(hapStatus))+(i-1)*length(hapStatus),strand.count[i]))
}
# expand the haplotype state dataframe according to the number of segments with different states
#repeat the rows based on the row indices defines in rep.rows
expanded.hapstates <- hapStrandStates[rep.rows]

# expand the probs datatable: repeat each row in the probs table #haps time
expanded.probs <- probs[sort(rep(1:nrow(probs),length(hapStatus)))]

#check if the strand states are the same in the two datatables
assert_that(all(expanded.probs$state == expanded.hapstates$state))

# combine the two datatables
new.probs <- cbind(expanded.probs, expanded.hapstates[,-"state"])

# compute dispersion parameters
new.probs <- add_dispPar(new.probs)

# compute NB haplotype likelihoods
new.probs[, nb_hap_ll := dnbinom(Wcn, size = disp_w, prob = nb_p)
                    *dnbinom(Ccn, size = disp_c, prob = nb_p)]

# not tested from here on
# averaging the nb probs of sister haplotypes
# , when haplotype specific strand states are not known

sister.hap.pos <- match(sister_haps, hapStatus)
sister.hap.pos[which(sister.hap.pos==1:length(hapStatus))]=NA
sister.hap.pos <- rep(sister.hap.pos, nrow(new.probs)/length(hapStatus))

if (!haplotypeMode)
{
  new_nb_hap_ll <- new.probs$nb_hap_ll
  new_nb_hap_ll[which(!is.na(sister.hap.pos))] <- (new.probs$nb_hap_ll
                      +new.probs$nb_hap_ll[sister.hap.pos])[which(!is.na(sister.hap.pos))]/2
  new.probs[, nb_hap_ll := new_nb_hap_ll]
}

# sorting the probs table based on segments
new.probs <- new.probs[order(sample, cell, chrom, from, to)]

# normalizing nb_hap_ll to 1 per sample, cell, and segment
new.probs[, nb_hap_ll := nb_hap_ll/sum(nb_hap_ll), by=.(sample, cell, chrom, from, to)]

# regularizing nb_hap_ll to set the min possible likelihood to a constant small number
# TODO complete this part

# add prior and compute the posteriori probs (add new columns)
# TODO complete this part

# normalizing again
new.probs[, nb_hap_ll := nb_hap_ll/sum(nb_hap_ll), by=.(sample, cell, chrom, from, to)]

# computing genotype likelihoods
nb_gt_ll <- new.probs$nb_hap_ll
nb_gt_ll[which(!is.na(sister.hap.pos))] <- (new.probs$nb_hap_ll
            +new.probs$nb_hap_ll[sister.hap.pos])[which(!is.na(sister.hap.pos))]
new.probs[, nb_gt_ll := nb_gt_ll]

# nb_gt_ll sums up to a very low number: double check

# converting to simple haplotype prob table

# calling SVs (It should be included in Sascha's code)