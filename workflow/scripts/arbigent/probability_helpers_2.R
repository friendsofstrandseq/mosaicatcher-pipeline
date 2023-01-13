#Results.

# First, load all the different bed files, put them together to a df
load_and_annotate_bedtable <- function(){
  
  # Get the different inversions
  basedir = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/results_overlord/used_beds/'
  invs_hom_link = paste(basedir,'HG00733_hom.bed', sep='')
  invs_het_link = paste(basedir,'HG00733_het.bed', sep='')
  invs_wt_link = paste(basedir,'HG00733_wt.bed', sep='')
  invs_shuf_link = paste(basedir,'iHGshuffsort.bed', sep='')
  
  invs_hom = read.table(invs_hom_link, header=F, sep='\t')
  invs_het = read.table(invs_het_link, header=F, sep='\t')
  invs_wt = read.table(invs_wt_link, header=F, sep='\t')
  invs_shuf = read.table(invs_shuf_link, header=F, sep='\t')
  
  # Annotate them
  invs_hom$gt = 'HOM inv'
  invs_het$gt = 'HET inv'
  invs_wt$gt = 'NO inv'
  invs_shuf$gt = 'Random segments'
  
  # Get them together in one df
  invs = rbind(invs_hom, rbind(invs_het, rbind(invs_wt, invs_shuf)))
  
  # Name columns
  colnames(invs) = c('chr','start','end','gt')
  return(invs)
}

make_call_consensus <- function(calls){
  
  library('dplyr')
  # just for demonstration 
  #calls = transform( calls, sv_call_name = sample(sv_call_name))
  
  calls_grouped = group_by(calls, chrom, start, end, sample )

  # Each group: count cells (should be the same), find the SV
  # with the highest count, calc af, and only keep the line that
  # contains this highest scoring SV. 
  top_sv_calls = count(calls_grouped, sv_call_name) %>%
    mutate(n_cells = sum(n)) %>%
    mutate(groupmax = max(n)) %>%
    mutate(af = n/n_cells) %>%
    top_n(1,af)
  
  return(top_sv_calls)
}

annotate_inv_with_callinfo <- function(invs, calls){
  for (row in 1:nrow(invs)){
    # If there is an SV call for this manual segment, give me the details for the first cell
    if (invs[row,]$start %in% calls$start){
      invs[row, c('sv_call_name', 'af')] = 
        calls[calls$start == invs[row,]$start,][1,c('sv_call_name', 'af')]
    } else {
      invs[row, c('sv_call_name', 'af')] =
        c('ref','1')
    }
  }
  return(invs)
}




# # Load SV calls
# call_file = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/results_overlord_fixed_thrice/a/simpleCalls_llr4_poppriorsTRUE_haplotagsTRUE_gtcutoff0_regfactor6.txt'
# calls_raw = read.table(call_file, sep='\t', header=1, stringsAsFactors = F)
# 
# # Get consensus call for each segment
# calls = make_call_consensus(calls_raw)
# 
# # Load inversion bedtable
# invs = load_and_annotate_bedtable()
# 
# # Annotate each manual segment with the consensus SV call
# invs = annotate_inv_with_callinfo(invs, calls)
# 
# invs$len = invs$end - invs$start
# invs$roundlen = round(log10(invs$len),1)
# 
# # make a plot #
# library(ggplot2)
# theme_set(theme_classic())
# 
# positions = c('HOM inv', 'HET inv', 'NO inv', 'Random segments')
# # Histogram on a Categorical variable
# g <- ggplot(invs, aes(gt)) +scale_fill_brewer(palette = "Spectral")
# 
# g + geom_bar(aes(fill=sv_call_name), width = 0.5) + 
#   theme(axis.text.x = element_text()) + 
#   scale_x_discrete(limits = positions) +
#   labs(title="SV calls", 
#        subtitle="") +
#   xlab('HG00733 Inversion status')
# 
# 
# # length related for the positive ones
# # make a plot #
# 
# invs_hom = invs[invs$gt == "HOM inv",]
# invs_het = invs[invs$gt == "HET inv",]
# invs_no = invs[invs$gt == "NO inv",]
# 
# library(ggplot2)
# theme_set(theme_classic())
# 
# #positions = c('HOM inv', 'HET inv', 'NO inv', 'Random segments')
# # Histogram on a Categorical variable
# g <- ggplot(invs_hom, aes(roundlen)) +scale_fill_brewer(palette = "Spectral")
# 
# g + geom_bar(aes(fill=sv_call_name), alpha=0.5) + 
#   labs(title="HOM inversion re-genotyping", 
#        subtitle="") +
#   xlab('log10 segment length')
# 
# #positions = c('HOM inv', 'HET inv', 'NO inv', 'Random segments')
# # Histogram on a Categorical variable
# g <- ggplot(invs_het, aes(roundlen)) +scale_fill_brewer(palette = "Accent")
# 
# g + geom_bar(aes(fill=sv_call_name), alpha=0.5, width=0.08) + 
#   labs(title="HET inversion re-genotyping", 
#        subtitle="") +
#   xlab('log10 segment length')
# 
# #positions = c('HOM inv', 'HET inv', 'NO inv', 'Random segments')
# # Histogram on a Categorical variable
# g <- ggplot(invs_no, aes(roundlen))# +scale_fill_brewer(palette = "Spectral")
# 
# g + geom_bar(aes(fill=sv_call_name), alpha=0.5, width=0.08) + 
#   labs(title="HET inversion re-genotyping", 
#        subtitle="") +
#   xlab('log10 segment length')
# 


