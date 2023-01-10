suppressMessages(library(data.table))
suppressMessages(library(assertthat))
suppressMessages(library(stringr))
source("workflow/scripts/arbigent_utils/mosaiclassifier_scripts/mosaiClassifier/mosaiClassifier.R")
source("workflow/scripts/arbigent_utils/mosaiclassifier_scripts/mosaiClassifier/haploAndGenoName.R")


#' Derives SV type with highest probability according in each cell and segment.
#' Output is a table with columns "sv_call_name" and "llr_to_ref"
#'
#' @author Sascha Meiers
#' @export
#'
makeSVCallSimple <- function(probs, llr_thr = 1, use.pop.priors = FALSE, use.haplotags = FALSE, genotype.cutoff = 0.0, bin.size, minFrac.used.bins = 0.5, manual.segs=FALSE) {

  assert_that(is.data.table(probs),
              "sample" %in% colnames(probs),
              "cell"   %in% colnames(probs),
              "chrom"  %in% colnames(probs),
              "start"  %in% colnames(probs),
              "end"    %in% colnames(probs),
              "haplo_name" %in% colnames(probs),
              "haplotype"  %in% colnames(probs),
              "nb_hap_ll"  %in% colnames(probs)) %>% invisible

  # make sure post-processing was done beforehands
  assert_that("nb_hap_pp" %in% colnames(probs)) %>% invisible

  # kick out the segments with a large fraction of blacklisted bins
  if(!manual.segs){
    probs <- probs[num_bins*bin.size/(end-start) >= minFrac.used.bins]
  }

  setkey(probs, chrom, start, end, sample, cell)

  if (use.pop.priors) {
    message('Applying population priors')
    probs[, pop.prior := sum(nb_hap_ll) , by = .(chrom, start, end, sample, haplotype)]
    probs[, pop.prior := pop.prior/sum(pop.prior) , by = .(chrom, start, end, sample, cell)]
    probs[, nb_hap_pp := nb_hap_pp*pop.prior]
  } else {
    message('Skipping population priors')
  }

  if (genotype.cutoff > 0.0) {
    message('Using genotype cutoff ', genotype.cutoff)
    probs[, nb_hap_pp_norm := nb_hap_pp/sum(nb_hap_pp) , by = .(chrom, start, end, sample, cell)]
    probs[, nb_hap_pp_norm_pop := sum(nb_hap_pp_norm) , by = .(chrom, start, end, sample, haplotype)]
    probs[, nb_hap_pp_norm_pop := nb_hap_pp_norm_pop/sum(nb_hap_pp_norm_pop) , by = .(chrom, start, end, sample, cell)]
    probs[, nb_hap_pp := ifelse(nb_hap_pp_norm_pop>genotype.cutoff, nb_hap_pp, 0.0)]
  } else {
    message('Skipping genotype cutoffs')
  }

  if (use.haplotags) {
    assert_that("haplotag.prob" %in% colnames(probs)) %>% invisible
    probs[!is.na(haplotag.prob), nb_hap_pp := nb_hap_pp*haplotag.prob ]
  }

  # annotate the ref_hom posterior probability per segment / cell
  probs[,
        ref_hom_pp := .SD[haplo_name == "ref_hom", nb_hap_pp[1]],
        by = .(chrom, start, end, sample, cell)]
        
  # order the different haplotype states based on their posterior prob. (nb_hap_pp)
  # and keep only the two most likely states
  probs[,
        rank := frank(-nb_hap_pp, ties.method = "random"), # random should never allow ties to get the same rank!
        by = .(chrom, start, end, sample, cell)]
  probs <- probs[rank <= 2]
  # Sort by rank within each group
  setkey(probs, chrom, start, end, sample, cell, rank)

  # select first and 2nd SV call
  probs[, `:=`(sv_call_name          = haplo_name[rank==1],
               sv_call_haplotype     = haplotype[rank==1],
               sv_call_name_2nd      = haplo_name[rank==2],
               sv_call_haplotype_2nd = haplotype[rank==2],
               llr_to_ref            = log(nb_hap_pp[rank==1]) - log(ref_hom_pp[rank==1]),
               llr_to_2nd            = log(nb_hap_pp[rank==1]) - log(nb_hap_pp[rank==2])),
       by = .(chrom, start, end, sample, cell)]
  probs <- probs[rank == 1]


  # Clean up table
  if(manual.segs){
    probs <- probs[, .(chrom, start, end, sample, cell, class, scalar, sv_call_name, sv_call_haplotype, sv_call_name_2nd, sv_call_haplotype_2nd, llr_to_ref, llr_to_2nd)]
  } else{
    probs <- probs[, .(chrom, start, end, sample, cell, class, scalar, num_bins, sv_call_name, sv_call_haplotype, sv_call_name_2nd, sv_call_haplotype_2nd, llr_to_ref, llr_to_2nd)]
  }

  return(probs[sv_call_name != "ref_hom" & llr_to_ref > llr_thr])
}

forceBiallelic <- function(probs, penalize_factor = 0.1)
{
  # Add reference probability as an extra column (same as in makeSimpleSVCalls)
  probs[,
        ref_hom_pp := .SD[haplo_name == "ref_hom", nb_hap_pp],
        by = .(chrom, start, end, sample, cell)]

  probs[, biall_hap_pp:=nb_hap_pp+ref_hom_pp]
  probs[haplo_name=="ref_hom", biall_hap_pp:=ref_hom_pp]
  setkey(probs,sample, chrom, start, end, haplotype)
  # computing aggregate biallelic probabilities
  probs[, agg_hap_pp := sum(log(biall_hap_pp)), by=.(sample, chrom, start, end, haplotype)]
  # adding the most Likely allele (other than reference) in the biallelic mode and creating a new column for that
  probs[, allele:=haplo_name[which.max(agg_hap_pp)], by=.(sample, chrom, start, end)]
  
  apply_prior <- function(probs, penalize_factor)
  {
    probs.new <- probs
    probs.new[haplo_name!=allele & haplo_name!="ref_hom", nb_hap_pp:=nb_hap_pp*penalize_factor]
    # normalization
    probs.new <- probs.new[, nb_hap_pp:=nb_hap_pp/(sum(nb_hap_pp)), by=.(sample, chrom, start, end, cell)]
    return(probs.new)
  }
  
  probs <- apply_prior(probs, penalize_factor)
}

# The following finction computes the most likely alternative allele for each segment only based on CC and WW cells
# and outputs the normalized prob table containing only the alt allele for each segment
#'
#' @author Maryam Ghareghani
#' @export
#'

getBiallelicLikelihoods <- function(probs, reg.factor)
{
  # remove complex SVs
  probs <- probs[geno_name != "complex"]

  # keep useful columns
  probs <- probs[, .(sample, cell, chrom, start, end, class, geno_name, nb_gt_ll)]
  probs <- unique(probs)

  # Add reference likelihoods as an extra column (same as in makeSimpleSVCalls)
  probs[,
        ref_hom_ll := .SD[geno_name == "ref_hom", nb_gt_ll],
        by = .(chrom, start, end, sample, cell)]

  probs[, biall_gt_ll:=nb_gt_ll+ref_hom_ll]
  probs[geno_name=="ref_hom", biall_gt_ll:=ref_hom_ll]
  setkey(probs, sample, chrom, start, end, geno_name)

  # get non-WC part of the probs table
  probs.not.wc <- probs[!class %in% c("WC", "CW")]
  # computing aggregate biallelic probabilities
  probs.not.wc[, agg_gt_ll := sum(log(biall_gt_ll)), by=.(sample, chrom, start, end, geno_name)]

  # adding the most Likely allele (other than reference) in the biallelic mode and creating a new column for that
  probs.not.wc[, alt_allele:=geno_name[which.max(agg_gt_ll)], by=.(sample, chrom, start, end)]
  
  # merge the new probs table (only for non-WC cells containing alt_allele) and the old probs table
  probs <- merge(probs, probs.not.wc, all.x=T, by=c("sample", "cell", "chrom", "start", "end", "geno_name", "class", "nb_gt_ll", "ref_hom_ll", "biall_gt_ll"))

  # remove agg_gt_ll column
  probs[, agg_gt_ll:=NULL]
  # define alt_allele for all cells such that each segment has a unique alt_allele over all cells
  probs[, alt_allele:=rep(alt_allele[!is.na(alt_allele)][1], .N), by=.(sample, chrom, start, end)]
  
  # keep only ref_hom and alt_allele for each segment
  probs <- probs[geno_name=="ref_hom" | geno_name==alt_allele]

  # normalize the prob table
  probs <- probs[, nb_gt_ll:=nb_gt_ll/sum(nb_gt_ll), by=.(sample, cell, chrom, start, end)]

  # remove segments with nan likelihood (TODO: look further into these cases)
  probs <- probs[!is.na(nb_gt_ll)]

  # keep only the alt_alle probs table
  probs <- probs[geno_name==alt_allele]

  # clean the probs table and remove unnessesary columns
  probs <- probs[, .(sample, cell, chrom, start, end, class, nb_gt_ll, alt_allele)]

  ### convert probs table to a wide matrix
  # define a column for SV event names
  probs[, event_name:=paste0(chrom, "_", start, "_", end, "_", alt_allele)]

  # regularize probs table
  # The underlying assumption is that there is a uniform distrubution (on ref and alt alleles) with reg.factor probability
  probs[, nb_gt_ll:=.((reg.factor/2)+nb_gt_ll*(1-reg.factor))]
  
  probs.mat <- dcast(probs, event_name~cell, value.var="nb_gt_ll")
  event_names <- probs.mat$event_name
  probs.mat <- as.matrix(probs.mat[,2:ncol(probs.mat)])
  rownames(probs.mat) <- event_names

  # set NA (uncertain) valuses in probs.mat matrix to 0.5
  probs.mat[which(is.na(probs.mat))] <- 0.5

  return(list(probs[, -"event_name"], probs.mat))
}

#' Derives SV type with highest likelihood in each cell and segment.
#' Output is a table with the new column "sv_call_name"
#'
#' @author Maryam Ghareghani
#' @export
#'
getMaxLikelihoodSV <- function(probs, haplotypeMode=F) {
  assert_that(is.data.table(probs),
              "sample" %in% colnames(probs),
              "cell"   %in% colnames(probs),
              "chrom"  %in% colnames(probs),
              "start"  %in% colnames(probs),
              "end"    %in% colnames(probs),
              "geno_name" %in% colnames(probs),
              "nb_gt_ll"  %in% colnames(probs)) %>% invisible

  if (haplotypeMode) {
      assert_that(is.data.table(probs),
              "haplo_name" %in% colnames(probs),
              "nb_hap_ll"  %in% colnames(probs)) %>% invisible
  }

  # remove complex SVs
  probs <- probs[geno_name!="complex"]

  # order different SVs based on their likelihoods
  if (haplotypeMode){
    setorder(probs, chrom, start, end, sample, cell, -nb_hap_ll)
  } else { 
    setorder(probs, chrom, start, end, sample, cell, -nb_gt_ll)
  }

  # keep only the most likely state
  probs <- probs[, head(.SD, 1), by=.(sample, cell, chrom, start, end)]

  if (haplotypeMode){
    colnames(probs)[which(colnames(probs)=="haplo_name")] <- "SV_class"
  } else { 
    colnames(probs)[which(colnames(probs)=="geno_name")] <- "SV_class"
  }
  
  return(probs[SV_class!="ref_hom"])
}

#' Computes CN likelihoods and Derives CN with highest probability according in each cell and segment.
#' Output is a table with columns "CN" and "CN_ll"
#'
#' @author Maryam Ghareghani
#' @export
#'
makeCNcall <- function(probs) {
  assert_that(is.data.table(probs),
              "sample" %in% colnames(probs),
              "cell"   %in% colnames(probs),
              "chrom"  %in% colnames(probs),
              "start"  %in% colnames(probs),
              "end"    %in% colnames(probs),
              "haplo_name" %in% colnames(probs),
              "haplotype"  %in% colnames(probs),
              "nb_hap_ll"  %in% colnames(probs)) %>% invisible

  # add CN column to probs table
  probs[, 
        CN:=haplo_code_to_geno_class(haplotype),
        by=haplotype]

  # compute CN likelihoods
  probs <- probs[,
                 .(CN_ll=sum(nb_hap_ll)),
                 by=.(sample,cell,chrom,start,end,CN)]

  # order the different copy numbers based on their likelihoods and keep only the most likely CN
  probs[,
        CN_rank := frank(-CN_ll, ties.method = "random"), # random should never allow ties to get the same rank!
        by = .(chrom, start, end, sample, cell)]
  probs <- probs[CN_rank < 2]
  probs[, CN_rank:=NULL]

  return(probs)
}
