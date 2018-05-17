library(data.table)
library(assertthat)
source("utils/mosaiClassifier/mosaiClassifier.R")

# delete: probs = readRDS("sv_probabilities/simulation1-50000/50000_fixed.medium/probabilities.Rdata")

makeSVCallSimple <- function(probs, llr_thr = 1) {

  # Do post-processing incl. priors + normalization + regularization
  probs = mosaiClassifierPostProcessing(probs)
  setkey(probs, chrom, start, end, sample, cell)

  # annotate the ref_hom posterior probability per segment / cell
  probs[,
        ref_hom_pp := .SD[haplo_name == "ref_hom", nb_hap_pp],
        by = .(chrom, start, end, sample, cell)]

  # order the different haplotype states based on their posterior prob. (nb_hap_pp)
  # and keep only the two most likely states
  probs[, rank := NULL]
  probs[,
        rank := frank(-nb_hap_pp, ties.method = "first"), 
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
               llr_to_2nd            = log(nb_hap_pp[rank==1]) - log(nb_hap_pp[rank==2]))]
  probs <- probs[rank == 1]


  # Clean up table
  probs <- probs[, .(chrom, start, end, sample, cell, class, sv_call_name, sv_call_haplotype, sv_call_name_2nd, sv_call_haplotype_2nd, llr_to_ref, llr_to_2nd)]

  return(probs[sv_call_name != "ref_hom" & llr_to_ref > llr_thr])
}


