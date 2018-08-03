suppressMessages(library(data.table))
suppressMessages(library(assertthat))
source("utils/mosaiClassifier/mosaiClassifier.R")


#' Derives SV type with highest probability according in each cell and segment.
#' Output is a table with columns "sv_call_name" and "llr_to_ref"
#'
#' @author Sascha Meiers
#' @export
#'
makeSVCallSimple <- function(probs, llr_thr = 1, use.pop.priors = FALSE, use.haplotags = FALSE, genotype.cutoff = 0.0, bin.size, minFrac.used.bins = 0.5) {

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
  probs <- probs[num_bins*bin.size/(end-start) >= minFrac.used.bins]

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
        ref_hom_pp := .SD[haplo_name == "ref_hom", nb_hap_pp],
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
  probs <- probs[, .(chrom, start, end, sample, cell, class, scalar, num_bins, sv_call_name, sv_call_haplotype, sv_call_name_2nd, sv_call_haplotype_2nd, llr_to_ref, llr_to_2nd)]

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
