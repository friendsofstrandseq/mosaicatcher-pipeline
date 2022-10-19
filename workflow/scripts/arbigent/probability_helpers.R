source("probability_helpers_2.R")
library(ggplot2)
library(ggbeeswarm)
library(reshape2)
library(dplyr)

# TODO NEEDS DESCRIPTION
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


#theme_set(theme_classic())
#' Return a tbl with median values for all haplotypes that were interested in. 
#' 
#' @param hapslist list of haplotypes (like c('0101','1010',...))
#' @param p_grouped  It's a tibble of probabilities derived from probabilites.Rdata, grouped by 
#' chrom, start, end

#' @author Wolfram Hoeps
#' @export
make_condensed_medianlist <- function(hapslist, p_grouped){
  counter = 1
  for (element in hapslist){
    if (counter == 1){
      cps = get_med_hap(element, p_grouped, mode='median')
    } else {
      cps = full_join(get_med_hap(element, p_grouped, mode='median'),cps, by=c("chrom","start","end"))
    }
    counter = counter + 1
  }
  return(cps)
}

#' Return a tbl with median values for all haplotypes that were interested in. 
#' 
#' @param hapslist list of haplotypes (like c('0101','1010',...))
#' @param p_grouped  It's a tibble of probabilities derived from probabilites.Rdata, grouped by 
#' chrom, start, end

#' @author Wolfram Hoeps
#' @export
make_condensed_sumlist <- function(hapslist, p_grouped){
  counter = 1
  for (element in hapslist){
    if (counter == 1){
      cps = get_med_hap(element, p_grouped, mode='sum')
    } else {
      cps = full_join(get_med_hap(element, p_grouped, mode='sum'),cps, by=c("chrom","start","end"))
    }
    counter = counter + 1
  }
  return(cps)
}

#' Return a tbl with median values for all haplotypes that were interested in. 
#' 
#' @param hapslist list of haplotypes (like c('0101','1010',...))
#' @param p_grouped  It's a tibble of probabilities derived from probabilites.Rdata, grouped by 
#' chrom, start, end

#' @author Wolfram Hoeps
#' @export
make_condensed_sumlist_probs <- function(hapslist, p_grouped){
  #haps_to_consider, pg_bulk_probs
  counter = 1
  for (element in hapslist){
    if (counter == 1){
      cps = get_med_hap(element, p_grouped, mode='prob')
    } else {
      cps = full_join(get_med_hap(element, p_grouped, mode='prob'),cps, by=c("chrom","start","end"))
    }
    counter = counter + 1
  }
  return(cps)
}

#' Search the grouped probabilites file for the entry of a certain haplotype (e.g. inv_h1) and
#' return the median across all cells. 
#' 
#' @param haplotype a string of the name of desired haplotype
#' @param p_grouped  It's a tibble of probabilities derived from probabilites.Rdata, grouped by 
#' chrom, start, end
#' @param haplo_or_geno
#' @param llr_limit cutoff for llr. default 10

#' @author Wolfram Hoeps
#' @param mode ['median', 'sum']. Today I think that sum is the best choice. 
#' @export
get_med_hap <- function(h_g_selection, pg, haplo_or_geno='haplotype', llr_limit=1e10, mode='median'){
  #if (!(h_g_selection %in% c('0101', '0202','0303','0404'))){
  sub_ps = pg[pg[[haplo_or_geno]]==h_g_selection,]
  #} else {
  #  sub_ps = pg[(pg[[haplo_or_geno]]==h_g_selection) & (pg$class %in% c('WW','CC')),]
  #}
  if (mode == 'median'){
    newcolname = paste('median_', h_g_selection, sep='')
    median_sub_ps = sub_ps %>% 
      mutate(new_mediancol = median(logllh)) %>%
      group_by(chrom, start, end) %>% 
      summarise(!!newcolname := mean(new_mediancol))
    median_sub_ps[[newcolname]][(median_sub_ps[[newcolname]]) > llr_limit] = llr_limit
    median_sub_ps[[newcolname]][(median_sub_ps[[newcolname]]) < -llr_limit/2] = -llr_limit/2
    
  } else if (mode == 'sum'){
    newcolname = paste('sum_', h_g_selection, sep='')
    # From this here, I removed a 'brob(dingnag)' from sum(exp(as.brob(logllh)))
    median_sub_ps = sub_ps %>% 
      mutate(new_mediancol = log(sum(exp(logllh)))) %>%
      group_by(chrom, start, end) %>% 
      summarise(!!newcolname := sum(new_mediancol))
    median_sub_ps[[newcolname]][(median_sub_ps[[newcolname]]) > llr_limit] = llr_limit
    median_sub_ps[[newcolname]][(median_sub_ps[[newcolname]]) < -llr_limit/2] = -llr_limit/2
  } else if (mode == 'prob'){
  newcolname = paste('sum_', h_g_selection, sep='')
  median_sub_ps = sub_ps %>% 
    mutate(new_mediancol = sum(realp)) %>%
    group_by(chrom, start, end) %>% 
    summarise(!!newcolname := mean(new_mediancol))
  #median_sub_ps[[newcolname]][(median_sub_ps[[newcolname]]) > llr_limit] = llr_limit
  #median_sub_ps[[newcolname]][(median_sub_ps[[newcolname]]) < -llr_limit/2] = -llr_limit/2
}
  return(median_sub_ps)
}

#' Remove from tbl entries that we do not like
#' I checked by hand the entries that systematically return LLR=0 in all cells.
# Then i filled the 'to_remove' list after checking the table by hand.
# pbad = pgi2[pgi2$logllh == 0,]
# unique_ones = (data.frame(pbad$haplotype, pbad$class))
# table(unique_ones)
#'
#' @param tb a grouped tbl
#' @author Someone from Stackexchange
#' @export
# to_remove = list(c('0011', c('WC', 'CW')),
#                  c('0101', c('WC', 'CW')),
#                   c('0020', c('WW', 'CC')),
#                   c('1100', c('WC', 'CW')),
#                   c('2000', c('WW', 'CC'))
#                  )
#filter_pg <- function(pgi_f){
#  pgi_filter = pgi_f %>% filter(
#    !((haplotype %in% c('0011', '0101','1100') & class %in% c('WC', 'CW')) |
#        (haplotype %in% c('0020', '2000') & class %in% c('WW', 'CC')))
#  )
#  
#  return(pgi_filter)
#}

# TODO: NEEDS DESCRIPTION
merge_inv_and_sdo <- function(invs_f, sd_overlap_f, bl_overlap_f){
  
  # Load overlap info, computer previously on bedtools
  sdo = read.table(sd_overlap_f, header=0, sep='\t'); colnames(sdo) = c('chr', 'start', 'end', 'sd_overlap')
  blo = read.table(bl_overlap_f, header=0, sep='\t'); colnames(blo) = c('chr', 'start', 'end', 'bl_overlap')
  
  # merge
  invs_m = full_join(invs_f, sdo)
  invs_m = full_join(invs_m, blo)
  
  # replace missing values with 0
  invs_m$sd_overlap[is.na(invs_m$sd_overlap)] = 0
  invs_m$bl_overlap[is.na(invs_m$bl_overlap)] = 0
  
  return(invs_m)
}

load_and_prep_pdf <- function(probs_file_f){
  pdf = data.frame(readRDS(probs_file_f))
  # Get LLR of each haplotype over reference
  # The implementation of this is not great, so be careful. 
  refs = pdf[pdf$haplotype=='1010',]$nb_hap_ll
  n_rep = dim(pdf)[1] / length(refs)
  pdf$nb_ref_ll = 0
  pdf$nb_ref_ll = rep(refs, each=n_rep) 
  pdf$logllh = log10(pdf$nb_hap_ll / pdf$nb_ref_ll)
  
  # pg 
  pg = group_by(pdf, chrom, start, end)
  pg = pg[!pg$class %in% '?',]
  
  # add len
  pg$len = pg$end - pg$start
  return(pg)
}


# TODO: NEEDS DESCRIPTION
give_me_inv <- function(probs_file_f){
  
  pg = load_and_prep_pdf(probs_file_f)
  
  # From the pg, remove observations that are meaningless. E.g. 0101 in WC
  #pg = filter_pg(pg)
  
  # Get a summary of that pg: median values per segment across cells
  haps_to_consider = unique(pg$haplotype)
  call_llhs = make_condensed_sumlist(haps_to_consider, pg)
  
  
  # Load inversion bedtable
  invs = load_and_annotate_bedtable()
  invs$len = invs$end - invs$start
  invs$roundlen = round(log10(invs$len),1)
  
  # Annotate each manual segment with the consensus SV call
  invs = full_join(invs, call_llhs)
  
  
  
  # For some reason an NA column can appear. Remove it. 
  invs = na.omit(invs)
  
  return(list(invs, pg, pdf))
}

# Can it be so easy??

#' Return the names of the haplotypes with the highest median LLR
#'
#' @param invs_f dataframe of inversions (_f for '_function')
#' @param pgi_f  pg stands for p_grouped. It's a tibble of probabilities 
#' derived from probabilites.Rdata (_f for '_function')
#' @param n   index of the inversion to be extracted
#' @author Wolfram Hoeps
#' @export
get_top_scoring_haplotypes <- function(invs_f, pgi_f, n, mode='sum'){
  relevant_inv = invs_f[invs_f$start == pgi_f$start[1],]
  if (mode=='median'){
    rel_inv_rel_cols = relevant_inv %>% select(starts_with('median_'))
    names_decreasing_withjunk = colnames(sort(rel_inv_rel_cols,decreasing = T))
    names_decreasing = substr(names_decreasing_withjunk,8,12)
    
    
  } else if (mode=='sum'){
    rel_inv_rel_cols = relevant_inv %>% select(starts_with('sum_'))
    names_decreasing_withjunk = colnames(sort(rel_inv_rel_cols,decreasing = T))
    names_decreasing = substr(names_decreasing_withjunk,5,9)
  }
  names_decreasing = names_decreasing[!names_decreasing %in% 'hetin']
  top_names = names_decreasing[0:n]
  return(top_names)
}

#' Return the names of the haplotypes with the highest median LLR
#'
#' @param call_llhs_f dataframe of inversions (_f for '_function')
#' @param pgi_f  pg stands for p_grouped. It's a tibble of probabilities 
#' derived from probabilites.Rdata (_f for '_function')
#' @param n   index of the inversion to be extracted
#' @author Wolfram Hoeps
#' @export
get_top_scoring_haplotypes_standalone <- function(call_llhs_f, pgi_f, n, mode='sum'){


  call_llhs_f2 = data.frame(call_llhs_f)
  relevant_inv = call_llhs_f2[(call_llhs_f2$start == pgi_f$start[1]) & (call_llhs_f2$end == pgi_f$end[1]),]
  if (mode=='median'){
    rel_inv_rel_cols = relevant_inv %>% select(starts_with('median_'))
    names_decreasing_withjunk = colnames(sort(rel_inv_rel_cols,decreasing = T))
    names_decreasing = substr(names_decreasing_withjunk,8,12)
    
    
  } else if (mode=='sum'){
    rel_inv_rel_cols = relevant_inv %>% select(starts_with('sum_'))
    names_decreasing_withjunk = colnames(sort(rel_inv_rel_cols,decreasing = T))
    if (rel_inv_rel_cols[names_decreasing_withjunk[1]] == 0){
      names_decreasing_withjunk = c('sum_1010', names_decreasing_withjunk[!names_decreasing_withjunk %in% c('sum_1010')])
    }
    names_decreasing = substr(names_decreasing_withjunk,5,9)
  }
  names_decreasing = names_decreasing[!names_decreasing %in% 'hetin']
  top_names = names_decreasing[0:n]
  return(top_names)
}


#' Extract a group from a grouped tbl
#'
#' @param x a grouped tbl
#' @param ngroup  index of the group to be extracted
#' @author Someone from Stackexchange
#' @export
select_group = function(x, ngroup) x %>%
  select(group_cols()) %>%
  distinct %>%
  ungroup %>%
  slice(ngroup) %>%
  { semi_join(x, .)}

translate = function(hapvec){
  t2 = hapvec
  for (i in 1:length(t2)){
    t2[i] = pgi[pgi$haplotype==hapvec[i],]$haplo_name[1]
  }
  return(t2)
}


maxN1 <- function(x, N=1){
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x,partial=len-N+1)[len-N+1]
}

maxN2 <- function(x, N=2){
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x,partial=len-N+1)[len-N+1]
}

attach_max_sec_name_to_call_llhs <- function(segs_llhs_f){ 

  segs_llhs = data.frame(segs_llhs_f)
  # if Hom or Het are tied for first place with something else, 
  # then we want to prefer Hom or Het. 
  segs_llhs$sum_0101 = segs_llhs$sum_0101# + 0.01
  segs_llhs$sum_0110 = segs_llhs$sum_0110# + 0.01
  segs_llhs$sum_0110 = segs_llhs$sum_0110# + 0.01
  segs_llhs[,'maxval'] = apply(segs_llhs %>% select(starts_with('sum_')), 1, maxN1)
  segs_llhs[,'secval'] = apply(segs_llhs %>% select(starts_with('sum_')), 1, maxN2)
  DF = (segs_llhs %>% select(starts_with('sum_')))
  segs_llhs[,'maxname'] = colnames(DF)[apply(DF,1,which.max)]
  segs_llhs[,'secname'] = colnames(DF)[apply(DF, 1, function(x)  which(x == sort(x, decreasing = TRUE)[2])[1])]
  segs_llhs[, "sum_hetinv_max"] <- apply(segs_llhs[, c('sum_0110', 'sum_1001')], 1, max)
  segs_llhs[, "interesting_minval"] <- apply(segs_llhs[, c('sum_hetinv_max','sum_0101','maxval')], 1, min)
  return(segs_llhs)
}



