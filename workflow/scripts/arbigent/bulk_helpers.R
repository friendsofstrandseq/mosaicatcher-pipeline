# Whoeps, 9th March 2021

#' @param pg_f
#' @author Wolfram Hoeps
#' @export 
#' 
calc_new_logllhs_singlecell <- function(pg_f){
  
  print('Re-calculating likelihoods after read count normalization. This will take up to several minutes.')
  # With new 'expected number of reads', recalculate dispersions W and C. Step 1
  
  pg_f1 = pg_f %>% group_by(chrom, start, end, class, haplotype) %>% mutate(  disp_w_prealpha=make_disp_w_prealpha(Wcn, expected, nb_p),
                                                                              disp_c_prealpha=make_disp_c_prealpha(Ccn, expected, nb_p))
  # With new 'expected number of reads', recalculate dispersions W and C. Step 2
  pg_f2 = pg_f1 %>% group_by(chrom, start, end, class, haplotype, Wcn, Ccn) %>% 
    mutate(disp_c_real = make_disp_c_real(Ccn, Wcn, disp_w_prealpha, disp_c_prealpha, expected, nb_p),
           disp_w_real = make_disp_w_real(Ccn, Wcn, disp_w_prealpha, disp_c_prealpha, expected, nb_p))
  
  # Round normalized count values 
  pg_f3 = pg_f2 %>% mutate(W = as.integer(round(W)),
                           C = as.integer(round(C)))
  
  # Now, calculate haplotype likelihoods
  pg_f4 = pg_f3 %>% mutate(nb_hap_ll = exp(recalc_nb_hap_ll(W, C, disp_c_real, disp_w_real, nb_p)))
  
  
  ### Process the nb_ref_llh ###
  # append the 'reference' likelihood to each line. For faster computation of likelihood ratio (next command)
  # Here a lazyer/more ugly part: We need to divide nb_hap_ll by reference so wg get logllh [i.e. logllh sv/ref].
  pg_f4$nb_ref_ll = 0
  for (cell in unique(pg_f4$cell)){
    refs = pg_f4$nb_hap_ll[(pg_f4$cell==cell) & (pg_f4$haplotype=='1010')]
    n_rep = dim(pg_f4[(pg_f4$cell==cell),])[1] / length(refs)
    pg_f4[(pg_f4$cell==cell),]$nb_ref_ll = rep(refs, each=n_rep) 
  }
  
  # log10 likelihood ratio vs reference
  pg_f4$logllh = log10(pg_f4$nb_hap_ll / pg_f4$nb_ref_ll)
  
  if (sum(pg_f4$logllh == 'Inf') > 0){
    pg_f4[pg_f4$logllh == 'Inf',]$logllh = 320
  }
  if (sum(pg_f4$logllh == 'Inf') > 0){
    pg_f4[pg_f4$logllh == 'Inf',]$logllh = 320
  }
  

  
  return(pg_f4)
  
  }

# Whoeps, July26th 2020

#' Function to turn a probabilities tibble into pseudo bulk. For this, we will
#' sum up the observed reads and expectations across cells of the same strand
#' state (WW,CC,WC,CW). Ultimately, this gives us a new probabilities tibble that
#' is structured like single cell, but contains only our 4 pseudo-bulk 'cells'.
#' 
#' Between the lines, we also have to recalculate a new NB model based on the 
#' Pseudo 'observations'. This is sone with the make_disp_... functions, and in the end 
#' with recalc_nb_hap_ll. The defintions for these is outsourced to other files though. 
#'
#' @param haps 
#' @param pg_f
#' @author Wolfram Hoeps
#' @export 
bulkify_pg <- function(haps, pg_f){

  #### NEEDS POLISHING AND DOCUMENTATION.
  
  # For the formula of nb_hap_ll, see mosaicatcher/utils/mosaiClassifier/getDispParAndSegType.R
  
  # We want to keep the stuff separated by start (i.e. segment), class (i.e. WW,WC,..) and haplytype (i.e. 1010, 1001, ...)
  # Here, count pseudobulk expected, W and C counts and make first step of disp_w and disp_c

  print('hi')
  pg2 = pg[ , !(names(pg) %in% c('logllh'))]
  # pre pre
  pg_f = na.omit(pg2)
  # Sometimes this omits only some haplotypes, while others survive. So let's say if one of the 70haps has na, we remove the whole inversion from the analysis.
  n = as.data.frame(rbind(table(pg_f$start)))
  bad_starts = colnames((n[as.numeric(n) %% 70 != 0]))
  pg_f = pg_f[!(pg_f$start %in% bad_starts),]
  
  
  
  pg_f = pg_f %>% group_by(class, haplotype) %>% mutate(Ccn = convert_hap(class, haplotype)[1], Wcn = convert_hap(class, haplotype)[2] )
  
  pg_f = merge_classes(pg_f)

  
  bulk_pre_1 = pg_f %>% group_by(chrom, start, end, class, haplotype) %>% mutate(expectedbulk = sum(expected), Wbulk = sum(W), Cbulk = sum(C),
                                                                     disp_w_prealpha=make_disp_w_prealpha(Wcn, expectedbulk, nb_p),
                                                                     disp_c_prealpha=make_disp_c_prealpha(Ccn, expectedbulk, nb_p))



  # Based on the class, make second step of disp_w, disp_c (...this is with this alpha stuff)
  # [W, 19th October]: removed the group_by here and added it to bulk_pre3 instead. Here, it caused an inexplicable error.

  bulk_pre_2 = bulk_pre_1 %>% #group_by(chrom, start, end, class, haplotype, Wcn, Ccn)  %>% 
                                      mutate(disp_c_real = make_disp_c_real(Ccn, Wcn, disp_w_prealpha, disp_c_prealpha, expectedbulk, nb_p),
                                        disp_w_real = make_disp_w_real(Ccn, Wcn, disp_w_prealpha, disp_c_prealpha, expectedbulk, nb_p))


  bulk_pre_2.5 = bulk_pre_2 %>% mutate(Wbulk = as.integer(round(Wbulk)),
                                       Cbulk = as.integer(round(Cbulk)))

  # Use our disp_w and disp_c, along with out pseudobulk expectation and counts, to calculate new hap_ll. 
  bulk_pre_3 = bulk_pre_2.5 %>% group_by(chrom, start, end, class, haplotype, Wcn, Ccn) %>% mutate(nb_hap_ll = recalc_nb_hap_ll(Wbulk, Cbulk, disp_c_real, disp_w_real, nb_p))

  # turn LLHs into probabilities, with fixed priors (1, 2 and 1/100, like in original mosaicatcher)
  bulk_pre_4 = bulk_pre_3 %>% mutate(nb_hap_prob = calc_prob(nb_hap_ll, haplotype),
                                     haplotype_prior = calc_prob(nb_hap_ll, haplotype, just_priors=T))

  # now remove all the unnecessary lines. we remain with pg_bulk with dimensions: [n invs x 2 cells x 70 haps]
  pg_bulk = bulk_pre_4 %>% group_by(chrom, start, end, haplotype, class, group) %>% slice(1)

  # Fix for small numbers. R does not allow numbers smaller than roughly 5e-324.
  # So it replaces them with 0. I.e. we have to replace the 0s with the smallest number then. 
  # leaving 4 magnitudes space :)
  #if (sum(pg_bulk$nb_hap_ll == 0) > 0){
  #  pg_bulk[pg_bulk$nb_hap_ll == 0,]$nb_hap_ll = 1e-320
  #  pg_bulk[pg_bulk$nb_hap_prob == 0,]$nb_hap_prob = 1e-320
  #}
  
  ### Process the nb_ref_llh ###
  # Here a lazyer/more ugly part: We need to divide nb_hap_ll by reference so wg get logllh [i.e. logllh sv/ref].
  pg_bulk$nb_ref_ll = 0
  for (cclass in c('CW','CC','WC','WW')){
    refs = pg_bulk$nb_hap_ll[(pg_bulk$class==cclass) & (pg_bulk$haplotype=='1010')]
    n_rep = dim(pg_bulk[(pg_bulk$class==cclass),])[1] / length(refs)
    pg_bulk[(pg_bulk$class==cclass),]$nb_ref_ll = rep(refs, each=n_rep) 
  }

  #pg_bulk$logllh1 = log10(pg_bulk$nb_hap_ll / pg_bulk$nb_ref_ll)
  pg_bulk$logllh = pg_bulk$nb_hap_ll - pg_bulk$nb_ref_ll
  
  if (sum(pg_bulk$logllh == 'Inf') > 0){
    pg_bulk[pg_bulk$logllh == 'Inf',]$logllh = 320
  }

  ### Process probabilities ###
  pg_bulk_probs = pg_bulk %>% group_by(chrom, start, end, haplotype) %>% mutate(product_probs = sum(log10(nb_hap_ll)) + log10(haplotype_prior))
  pg_bulk_probs = pg_bulk_probs %>% group_by(chrom, start, end, haplotype) %>% slice(1) %>% group_by(chrom, start, end) %>% mutate(summ = sum(10**product_probs), realp = (product_probs)-log10(summ))
  #pg_bulk = pg_bulk %>% group_by(start, end, class) %>% mutate(prob = nb_hap_prob/sum(nb_hap_prob))

  # Rename our newly created frankensteinbulkcell
  pg_bulk$cell = paste0('bulkcell_', pg_bulk$class)
  pg_bulk_probs$cell = paste0('bulkcell_', pg_bulk_probs$class)
  

  return(list(pg_bulk, pg_bulk_probs))
  
}



merge_classes <- function(pg_f2){
  
  # We want to count WW+CC and WC+CW cells together. So there is some pre-preprocessing needed.
  # from https://stackoverflow.com/questions/7746567/how-to-swap-values-between-two-columns
  # pre_pre 1) counts
  pg_f2 <- transform(pg_f2, W = ifelse(class == 'WW', C, W), C = ifelse(class == 'WW', W, C))
  pg_f2 <- transform(pg_f2, W = ifelse(class == 'WC', C, W), C = ifelse(class == 'WC', W, C))
  
  # pre_pre 2) copy numbers
  pg_f2 <- transform(pg_f2, Wcn = ifelse(class == 'WW', Ccn, Wcn), Ccn = ifelse(class == 'WW', Wcn, Ccn))
  pg_f2 <- transform(pg_f2, Wcn = ifelse(class == 'WC', Ccn, Wcn), Ccn = ifelse(class == 'WC', Wcn, Ccn))

  # pre_pre 3) class names
  pg_f2$class[pg_f2$class=='WW'] = 'CC'
  pg_f2$class[pg_f2$class=='WC'] = 'CW'
  
  return(pg_f2)
}

recalc_nb_hap_ll <- function(Wbulk, Cbulk, disp_c_real, disp_w_real, nb_p){

  # Use all our parameters to calc nb_hap_ll. 
  # Again, see mosaicatcher/utils/mosaiClassifier/getDispParAndSegType.R and 
  # mosaicatcher/utils/mosaiClassifier/mosaiClassifier.R to understand. 
  #nb_hap_ll = exp(dnbinom(Wbulk, size=disp_w_real, prob=nb_p, log=T) + dnbinom(Cbulk, size=disp_c_real, prob=nb_p, log=T))
  nb_hap_ll = (dnbinom(Wbulk, size=disp_w_real, prob=nb_p, log=T) + dnbinom(Cbulk, size=disp_c_real, prob=nb_p, log=T))
  
  return(nb_hap_ll)
}


calc_prob <- function(nb_hap_ll, haplotype,just_priors=F){
  # We are using bayes: p(SV|data) = (p(data|SV) * p(SV)) / p(data).
  # we don't know p(data), but it is just a normalization factor 
  # so we don't need it. 
  # p(data|SV) is simply nb_hap_ll. And p(SV) we can choose freely (that's our prior!).
  # We choose the prior to be 1 for invs, dups and dels, 2 for reference (to be conservative)
  # and 0.01 for anything else (aka complex rearrangements)
  haps = unique(haplotype)
  highpriorhaps = c('0000','0010','1000','0101','0110',
                    '1001','1011','1110','2020','1020',
                    '2010','1111')
  highestpriorhaps = '1010'
  priortable = data.frame(probs = rep(1/100,length(haps)))
  row.names(priortable) = haps
  priortable[highpriorhaps,] = 1
  priortable[highestpriorhaps,] = 2
  if (just_priors){
    return(priortable[haplotype,])
  }
  nb_hap_prob = nb_hap_ll * priortable[haplotype,]
  return(nb_hap_prob)
}

make_disp_w_real <- function(Ccn, Wcn, disp_w_prealpha, disp_c_prealpha, expected, nb_p, alpha=0.05, scalar=1 ){
  # Take second step of calculating disp_w
  # Again, see mosaicatcher/utils/mosaiClassifier/getDispParAndSegType.R and 
  # mosaicatcher/utils/mosaiClassifier/mosaiClassifier.R to understand. 
  
  ### Wcn and Ccn are zero? ###
  if (all(Ccn==0) & all(Wcn==0)){
    disp_w = scalar * expected * (nb_p / (1-nb_p)) * alpha * 0.5
  } else if (all(Ccn==0) & all(Wcn>0)) {
    disp_w = disp_w_prealpha * (1-alpha)
  } else if (all(Ccn>0) & all(Wcn==0)) {
    disp_w = disp_c_prealpha * alpha
  } else {
    disp_w = disp_w_prealpha
  }
  return(as.numeric(disp_w))
}

make_disp_c_real <- function(Ccn, Wcn, disp_w_prealpha, disp_c_prealpha, expected, nb_p, alpha=0.05, scalar=1 ){
  # Take second step of calculating disp_c
  # Again, see mosaicatcher/utils/mosaiClassifier/getDispParAndSegType.R and 
  # mosaicatcher/utils/mosaiClassifier/mosaiClassifier.R to understand. 
  # print(Ccn)
  # print(Wcn)
  # print("go")
  ### Wcn and Ccn are zero? ###
  if (all(Ccn==0) & all(Wcn==0)){
    disp_c = scalar * expected * (nb_p / (1-nb_p)) * alpha * 0.5
  } else if (all(Ccn==0) & all(Wcn>0)) {
    disp_c = disp_w_prealpha * (alpha)
  } else if (all(Ccn>0) & all(Wcn==0)) {
    disp_c = disp_c_prealpha * (1-alpha)
  } else {
    disp_c = disp_c_prealpha
  }  
  return(as.numeric(disp_c))
}


make_disp_w_prealpha <- function(Wcn, expected, nb_p, alpha=0.05, scalar=1){
  # Take first step of calculating disp_w
  # Again, see mosaicatcher/utils/mosaiClassifier/getDispParAndSegType.R and 
  # mosaicatcher/utils/mosaiClassifier/mosaiClassifier.R to understand. 
  disp_w_prealpha = scalar * expected * nb_p / (1-nb_p) * Wcn * 0.5
  return(disp_w_prealpha)
}

make_disp_c_prealpha <- function(Ccn, expected, nb_p, alpha=0.05, scalar=1){
  # Take first step of calculating disp_c
  # Again, see mosaicatcher/utils/mosaiClassifier/getDispParAndSegType.R and 
  # mosaicatcher/utils/mosaiClassifier/mosaiClassifier.R to understand. 
  disp_w_prealpha = scalar * expected * nb_p / (1-nb_p) * Ccn * 0.5
  return(disp_w_prealpha)
}



convert_hap <- function(class, hap_f){
  # Ugly but it runs. 
  # 0=Crick
  # 1=Watson
  strands = strsplit(class,'')[[1]]
  strands[strands=='C'] = 0
  strands[strands=='W'] = 1
  strands=as.numeric(strands)
  haps = c(rep(strands[1], as.numeric(strsplit(hap_f,'')[[1]])[1]),
           rep(1-strands[1], as.numeric(strsplit(hap_f,'')[[1]])[2]),
           rep(strands[2], as.numeric(strsplit(hap_f,'')[[1]])[3]),
           rep(1-strands[2], as.numeric(strsplit(hap_f,'')[[1]])[4]))
  haps[haps==0] = 'C'
  haps[haps==1] = 'W'
  Ccn = length(grep('C', haps))
  Wcn = length(grep('W', haps))
  #return(Ccn)
  return(c(Ccn,Wcn))
}

