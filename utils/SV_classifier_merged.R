# parameters:
maximumCN <- 4
haplotypeMode <- F
regularizationFactor <- 1e-10

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
  probs <- prob.tab
  probs[, disp_w := scalar*nb_r*Wcn*(to-from+1)*(0.5)]
  probs[, disp_c := scalar*nb_r*Ccn*(to-from+1)*(0.5)]
  
  # setting both W and C dispersion params for CN0 to alpha
  cn0_ridx <- which(probs$disp_c==0 & probs$disp_w==0)
  probs[cn0_ridx, disp_w:=scalar*nb_r*alpha*(to-from+1)]
  probs[cn0_ridx, disp_c:=scalar*nb_r*alpha*(to-from+1)]
  
  # rescaling W and C dispersion params for Wcn=0 cases with the parameter alpha
  Wcn0_ridx <- which(probs$disp_w==0)
  probs[Wcn0_ridx, disp_w:=disp_c*alpha*2]
  probs[Wcn0_ridx, disp_c:=disp_c*(1-alpha)*2]
  
  # rescaling W and C dispersion params for Ccn=0 cases with the parameter alpha
  Ccn0_ridx <- which(probs$disp_c==0)
  probs[Ccn0_ridx, disp_c:=disp_w*alpha*2]
  probs[Ccn0_ridx, disp_w:=disp_w*(1-alpha)*2]
  
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
    
    segType = paste0(stringr::str_dup("W",Nw),stringr::str_dup("C",Nc))
  }
  
  #list(segType, Wcn = Nw, Ccn = Nc)
  c(Nw, Nc)
}

#' returns the hapotyoe name
#' 
#' @param hap.code The haplotype coding
#' @author Maryam Ghareghani
#' @export
#' 

get_hap_name <- function(hap.code)
{
  hap.codes <- c("1010", "0001", "1000", "0000", "0110", "1001", "0101", "2010", "1020", "2020")
  hap.names <- c("ref_hom", "del_h1", "del_h2", "del_hom", # ref and del
                                  "inv_h1", "inv_h2", "inv_hom", # inv
                                  "dup_h1", "dup_h2", "dup_hom") # dup
  
  hap.idx <- match(hap.code, hap.codes)
  
  if (!is.na(hap.idx))
  {
    return(hap.names[hap.idx])
  }
  
  return("complex")
}

#' converts the haplotype to the genotype name
#' 
#' @param hap.name The haplotype name
#' @author Maryam Ghareghani
#' @export
#' 

haplo_to_geno_name <- function(hap.name)
{
  geno.name <- gsub("h1","het",hap.name)
  geno.name <- gsub("h2","het",geno.name)
  
  return(geno.name)
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
# naming third and forth columns
colnames(hapStrandStates)[3:4]=c("Wcn", "Ccn")
# adding haplotype and genotype name columns
hapStrandStates[,haplo_name:=.(sapply(haplotype, get_hap_name))]
hapStrandStates[,geno_name:=.(sapply(haplo_name, haplo_to_geno_name))]
# sort based on state
setkey(hapStrandStates, state)
##### COMMENT: this part probably can be done in a better way using datatable
# I don't know how to merge these two datatables (probs and hapStrandStates)
# took the only solution that I have in mind

# kick out the segs with sces
#probs <- probs[state!="sce"]

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

# expand the probs datatable: First sort the rows based on states, then repeat each row in the probs table #haps time
expanded.probs <- probs[order(state)][sort(rep(1:nrow(probs),length(hapStatus)))]

#check if the strand states are the same in the two datatables
assert_that(all(expanded.probs$state == expanded.hapstates$state))

# combine the two datatables
probs <- cbind(expanded.probs, expanded.hapstates[,-"state"])
###########
# compute dispersion parameters
probs <- add_dispPar(probs)

# compute NB haplotype likelihoods
probs[, nb_hap_ll := dnbinom(Wcn, size = disp_w, prob = nb_p)
                    *dnbinom(Ccn, size = disp_c, prob = nb_p)]

# computing sister haplotype (haplotype with the same genotype) for each haplotype
sister.haps <- sapply(hapStatus, sisterHaplotype)
sister.hap.pos <- match(sister.haps, hapStatus)
# compute the set of symmetric haplotypes (haplotypes that are equal to their sister haplotype)
symmetric.haps <- hapStatus[which(sister.haps==hapStatus)]
# testing:
# test <- probs[haplotype=="1000"|haplotype=="0010",]
# test[, diff_ll:=.SD$nb_hap_ll[2]-.SD$nb_hap_ll[1],by=.(sample, chrom, cell, from, to)]
# unique(test$diff_ll)
# It did not work for the tested data, all pairs of sister haplotypes turned to have the same ll
# I should test the following part later

# averaging the nb probs of sister haplotypes, when haplotype specific strand states are not known
# adding genotype likelihoods, if haplotype mode is false
if (!haplotypeMode)
{
  probs[,nb_hap_ll:=.((nb_hap_ll+nb_hap_ll[sister.hap.pos])/2), by=.(sample, cell, chrom, from, to)]
  
  # computing genotype likelihoods
  probs[,nb_gt_ll:=.(nb_hap_ll+nb_hap_ll[sister.hap.pos]), by=.(sample, cell, chrom, from, to)]
  # deviding the gt likelihoods of symmetric haplotypes by 2
  probs[haplotype %in% symmetric.haps, nb_gt_ll:=.(nb_gt_ll/2)]
}

# TODO export the prob table to some output file

### Post processing

# testing if there are some segments with zero probability for all haplotypes
segs_max_hap_nb_probs <- probs[, .(sample, chrom, cell, from, to, max_nb_hap_ll=rep(max(nb_hap_ll), .N)), 
                              by=.(sample, chrom, cell, from, to)]
# there are some 0s (0.3% of segments in the tested example) ???(how to treat them)
# segs_max_hap_nb_probs[max_nb_hap_ll==0] # non-empty

# add prior probs to the table
# TODO set the class of the prior column to numeric, it set the prior of complex haplos to z now
probs[,prior:=1.0L]
probs[haplo_name=="ref_hom",prior:=2L]
probs[haplo_name=="complex",prior:=.(prior/100)]

# compute the posteriori probs (add new columns)
probs[,nb_hap_pp:=.(nb_hap_ll*prior)][,nb_gt_pp:=.(nb_gt_ll*prior)]

# set a uniform prob on sce segs and the segs_max_hap_nb_probs=0
# TODO test this part
probs[segs_max_hap_nb_probs$max_nb_hap_ll==0,nb_hap_pp:=1L]
probs[state=="sce", nb_hap_pp:=1L]

# normalizing nb_hap_pp to 1 per sample, cell, and segment
probs[, nb_hap_pp := nb_hap_pp/sum(nb_hap_pp), by=.(sample, cell, chrom, from, to)]
# some NANs are produced here (It will be resolved after prev #TODO)

# regularizing nb_hap_ll to set the min possible likelihood to a constant small number
# TODO test this part
probs[,nb_hap_pp:=.((regularizationFactor/length(hapStatus))+nb_hap_pp*(1-regularizationFactor))]

# computing genotype posterioris
probs[,nb_gt_pp:=.(nb_hap_pp+nb_hap_pp[sister.hap.pos]), by=.(sample, cell, chrom, from, to)]
# deviding the gt likelihoods of symmetric haplotypes by 2
probs[haplotype %in% symmetric.haps, nb_gt_pp:=.(nb_gt_ll/2)]

# converting to simple haplotype prob table
simp.probs <- probs[haplo_name!="complex"]

# normalizing hap and gt probs in simp.probs
probs[, nb_hap_pp := nb_hap_pp/sum(nb_hap_pp), by=.(sample, cell, chrom, from, to)]
probs[, nb_gt_pp := nb_hap_pp/sum(nb_gt_pp), by=.(sample, cell, chrom, from, to)]

# dcasting: converting the table from long to wide format based on the haplotypes

# calling SVs (It should be included in Sascha's code)