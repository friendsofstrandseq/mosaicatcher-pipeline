# Whoeps, 07th Jan 2021
# Making a large overview over the results from the arbigent folder. 
# I'm giving myself 2h to make this nice today. 

# Input: callmatrix from clean_genotype.R
# Input: a csv from david from which to extract samplenames
# Output: a matrix with added entries:
#   Filter - Pass, NoReadsPass, MendelFail, FalsePositive, lowconf, AlwaysComplex
#   Mapability - percent?
#   nhom, nhet, nref, nnoreads, ncomplex
#   mendel 1/0

# Load libraries
library(stringr)
library(dplyr)
library(pheatmap)
library(matrixStats)
library(reshape2)
# Define functions
extract_used_samples <- function(csvlink){
  # Load david's genotypes
  
  dgt = read.table(david_list, header=T, sep=',', stringsAsFactors = F)
  
  # define entries to keep
  cols_to_keep = c('genoT')
  
  # filter
  dgtf = dgt %>%   select(matches(paste(cols_to_keep, collapse="|")))
  
  # rename columns
  colnames(dgtf) = str_replace(colnames(dgtf),'genoT_','')
  
  return(colnames(dgtf))
  
}

cut_cm_to_samples <- function(cm_f, used_samples_uniq_f, n_overjump, include_children){
  if (include_children){
    children_ids = c('HG00514','HG00733','NA19240')
  } else {
    children_ids = NULL
  }
  # Change callmatrix samplenames
  for (i in length(colnames(cm_f)):n_overjump){
    if (colnames(cm_f)[i]=='HG002'){
      colnames(cm_f)[i] = 'NA24385'
      next
    } else if (colnames(cm_f)[i] %in% children_ids){
      next
    }
    newlabel = used_samples_uniq_f[grepl(substr(colnames(cm_f)[i],3,8),used_samples_uniq_f)]
    print(newlabel)
    if (length(newlabel)>0){
      colnames(cm_f)[i] = newlabel
    } else {
      #print(paste0('Dropped sample ',colnames(cm)[i]))
      cm_f[,colnames(cm_f)[i]] = NULL
    }
    newlabel = NULL
  }

  return(cm_f)
}

count_homhetrefetc <- function(cm_f){
  cm_f$nhom = rowSums2(cm_f=='1|1') + rowSums2(cm_f=='1|1_lowconf')
  cm_f$nhet = rowSums2(cm_f=='1|0') + rowSums2(cm_f=='0|1') + rowSums2(cm_f=='1|0_lowconf') + rowSums2(cm_f=='0|1_lowconf')
  cm_f$nref = rowSums2(cm_f=='0|0') + rowSums2(cm_f== '0|0_lowconf')
  cm_f$nnoreads = rowSums2(cm_f=='noreads')
  cm_f$ncomplex = length(used_samples_uniq) - (cm_f$nhom + cm_f$nhet + cm_f$nref + cm_f$nnoreads)
  return(cm_f)
}

apply_filter <- function(cm_f, used_samples_uniq){
  # Judgement hour: filter
  cm_f$gt_events = cm_f$nhom + cm_f$nhet + cm_f$ncomplex 
  
  # Criteria hour
  # If we have < 5 bins, we do not attempt to filter
  
  # Failure modes
  # FalsePositive: no nonref events confirmed
  # AlwaysComplex: no noncomplex events confirmed
  
  # Explicit pass modes
  # Undercalled: CR is below 0.5, and we observe more events in GT than in the discovery list.
  # Pass: Events match perfectly 
  # NoReadsPass: Events match perfectly, difference can be explained by nnoreads (?)
  # PassMulticallers 
  
  # Judgement hour
  cm_f$verdict = 'UNK'
  bincutoff = 5
  cm_f[(cm_f$valid_bins <= bincutoff),]$verdict = 'lowconf'
  cm_f[(cm_f$gt_events == 0) & (cm_f$valid_bins > bincutoff),]$verdict = 'NoEvents'
  cm_f[(cm_f$mendelfails > 0) & (cm_f$valid_bins > bincutoff),]$verdict = 'MendelFail'
  cm_f[((cm_f$nhom + cm_f$nhet + cm_f$nref) == 0) & (cm_f$valid_bins > bincutoff),]$verdict = 'AlwaysComplex'
  
  cm_f[(cm_f$nnoreads == length(used_samples_uniq)),]$verdict = 'NoReadsPass'
  #cm_f[(cm_f$gt_events == cm_f$eventsclaimed),]$verdict = 'Pass'
  cm_f[cm_f$verdict == 'UNK',]$verdict = 'Pass'
  #cm_f[cm_f$CALLERSET_LIST %in% c('PAV,SSEQAUTO','PAV,SSEQAUTO,BIONANO','SSEQAUTO,BIONANO'),]$verdict = 'Pass_multicallers'
  return(cm_f)
}

#' Make a vector to check inheritance plaisibilities
#' @return child_expect, a vector describing excpected child genotypes given the parents.
#' @author Wolfram Hoeps
make_child_expect_vector <- function(){
  # ok who cares let's do it the hard way.
  #parents = paste(p1, p2)
  child_expect <- vector(mode="list")
  child_expect[['0|0 0|0']] = c('0|0')
  child_expect[['0|0 0|1']] = c('0|0','0|1')
  child_expect[['0|0 1|1']] = c('0|1')
  child_expect[['0|1 0|0']] = c('0|0','0|1')
  child_expect[['0|1 0|1']] = c('0|0','0|1','1|1')
  child_expect[['0|1 1|1']] = c('0|1','1|1')
  child_expect[['1|1 0|0']] = c('0|1')
  child_expect[['1|1 0|1']] = c('0|1','1|1')
  child_expect[['1|1 1|1']] = c('1|1')
  return (child_expect)
}

#' Test a gt for validity
#' 
test_mendel <- function(ce, gt_parent1, gt_parent2, gt_child){
  gt_parent1 = substr(gt_parent1,1,3)
  gt_parent2 = substr(gt_parent2,1,3)
  gt_child = substr(gt_child, 1,3)
  gt_parent1 = gsub('1\\|0','0\\|1', gt_parent1)
  gt_parent2 = gsub('1\\|0','0\\|1', gt_parent2)
  gt_child   = gsub('1\\|0','0\\|1', gt_child)
  
  valid_gts = c('0|0','0|1','1|1')
  c1 = gt_parent1 %in% valid_gts
  c2 = gt_parent2 %in% valid_gts
  c3 = gt_child %in% valid_gts
  #print(c1)
  #print(c2)
  #print(c3)
  #return(gt_parent2)
  if (c1){
    if (c2){
      if (c3){
        valid = gt_child %in% ce[[paste(gt_parent1, gt_parent2)]]
        #print(valid)
        return(valid)   
      }
    }
  }
  return('CPX')
}


add_mendelfails <- function(cm_f){
  ce = make_child_expect_vector()
  callmatrix = cm_f
  callmatrix$mendel1 = 'UNK'
  callmatrix$mendel2 = 'UNK'
  callmatrix$mendel3 = 'UNK'
  callmatrix$mendel4 = 'UNK'
  
  for (row in 1:nrow(callmatrix)){
    #callmatrix[row,]$mendel = 
    callmatrix[row,]$mendel1 = (test_mendel(ce, callmatrix[row,]$NA19238, callmatrix[row,]$NA19239,callmatrix[row,]$NA19240 ))
    callmatrix[row,]$mendel2 = (test_mendel(ce, callmatrix[row,]$HG00512, callmatrix[row,]$HG00513,callmatrix[row,]$HG00514 ))
    callmatrix[row,]$mendel3 = (test_mendel(ce, callmatrix[row,]$HG00731, callmatrix[row,]$HG00732,callmatrix[row,]$HG00733 ))
    #callmatrix[row,]$mendel4 = as.logical(test_mendel(ce, callmatrix[row,]$GM19650, callmatrix[row,]$HG00864,callmatrix[row,]$HG03371 ))
    
  }
  
  callmatrix$mendelfails = rowSums(callmatrix[,c('mendel1','mendel2','mendel3')] == 'FALSE')
  return(callmatrix)
}


# parameters
run_unsupervised = T
run_comparison = F
include_children = T

# inputs
david_list = '~/PhD/projects/huminvs/mosaicatcher/bed_factory/revision/david_n35_323/nonred_inversions_n35_genotypes.csv'
used_samples_uniq = extract_used_samples(david_list)

# hardcoded parameters
n_overjump = 7 # This you have to check in cm. How many cols at the beginning are not samples? (chrom, start, end, id,...)

if (run_unsupervised){
  ##########################################################################
  ### THIS IS BASED ON THE GT RESULTS ALONE. YOU COULD SAY UNSUPERVISED. ###
  ##########################################################################
  
  # Get callmatrix
  cm = callmatrix
  
  # cut down cm to the samples we are interested in
  cm = cut_cm_to_samples(cm, used_samples_uniq, n_overjump, include_children)
  
  # Factor char stuff
  cm[] <- lapply(cm, as.character)
  
  # stratify entries with 0 valid bins. 
  cm[cm$valid_bins==0,used_samples_uniq] = 'noreads'
  
  # Count hom, het, ref, noreads and complex
  cm = count_homhetrefetc(cm)
  
  # Calc mapability
  cm$mapability = (as.numeric(cm$valid_bins)/10) / as.numeric(cm$len)
  cm$valid_bins = as.numeric(cm$valid_bins)
  
  # Mendel
  cm = add_mendelfails(cm)
  
  # Filter
  cm = apply_filter(cm, used_samples_uniq)
  cm[] <- lapply(cm, as.character)
}


#############################################################################################
### THIS IS A COMPARISON BETWEEN REGENOTYPER AND ANOTHER CALLSE (IN THIS CASE FROM DAVID) ###
#############################################################################################

load_davidtable <- function(david_list){
  dgt = read.table(david_list, header=T, sep=',', stringsAsFactors = F)
  
  # define entries to keep
  cols_to_keep = c('seqnames','start','end','width','genoT')
  
  # filter
  dgtf = dgt %>%   select(matches(paste(cols_to_keep, collapse="|")))
  
  # rename columns
  colnames(dgtf) = str_replace(colnames(dgtf),'genoT_','')
  
  # we have to do some more renamings.
  colnames(dgtf) = str_replace(colnames(dgtf),'A$','')
  colnames(dgtf) = str_replace(colnames(dgtf),'B$','')
  
  return(dgtf)
}

append_homhetetc <- function(cm_f){
    cm_f$nhom_d = rowSums2(cm_f=='HOM')
    cm_f$nhet_d = rowSums2(cm_f=='HET')
    cm_f$nref_d = rowSums2(cm_f=='REF') 
    cm_f$nnoreads_d = rowSums2(cm_f=='lowReads')
    return(cm_f)
}
  
rename_gter_gts <- function(df){
  df$gtyper_simpler = 'complex'
  df[df$genotyper %in% c('1|1', '1|1_lowconf'),]$gtyper_simpler = 'HOM'
  df[df$genotyper %in% c('1|0', '1|0_lowconf','0|1', '0|1_lowconf'),]$gtyper_simpler = 'HET'
  df[df$genotyper %in% c('0|0', '0|0_lowconf'),]$gtyper_simpler = 'REF'
  df[as.numeric(df$valid_bins) < 5,]$gtyper_simpler = 'Lowconf'
  df[df$genotyper %in% c('noreads'),]$gtyper_simpler = 'zeroreads'
  return(df)
}

# load davidtable
dgtf = load_davidtable(david_list)
dgtf_numbers = append_homhetetc(dgtf)
dgtf_numbers <- as.data.frame(lapply(dgtf_numbers, as.character))

# Get the regenotyper data, adjust column names
cmm2 = cut_cm_to_samples(callmatrix, used_samples_uniq, 7, include_children)
mygts = melt(data.frame(cmm2),  id.vars = c("chrom", "start",'end','ID','len','valid_bins'))
colnames(mygts) = str_replace(colnames(mygts),'chrom','seqnames')
colnames(mygts) = str_replace(colnames(mygts),'variable','sample')

# melt dgtf
davidgts = melt(dgtf, id.vars = c('seqnames','start','end','width'))
colnames(davidgts) = str_replace(colnames(davidgts),'variable','sample')

# join things
joined = inner_join(mygts, davidgts, by=c('seqnames','start','end','sample'))
colnames(joined) = c('seqnames','start','end','ID','width','valid_bins','sample','genotyper','ignore','david')

# rename stuff again
joined = rename_gter_gts(joined)

# Count matches
joined$match = (joined$gtyper_simpler == joined$david)

# Cast into matrix form
jmat = cast(joined, seqnames+start+end+width~sample, value='match')

# Count matches per inversion, and in pct.
jmat$match = rowSums(jmat==T)
jmat$nomatch = rowSums(jmat==F)
jmat$matchpct = round((jmat$match/(jmat$match + jmat$nomatch)), 3)
jmat[] <- lapply(jmat, as.character)
cm2 = left_join(cm, jmat[,c('seqnames', 'start', 'end', 'match', 'nomatch', 'matchpct')])
cm3 = left_join(cm2, dgtf_numbers[,c('seqnames','start','end','nhom_d','nhet_d','nref_d','nnoreads_d')])
cm4 = cm3[!names(cm3) %in% used_samples_uniq]

# Sort into better order
cm5 = cm4[,c('verdict','chrom','start','end','ID','len','valid_bins','mapability','mendelfails','nhom','nhet','nref','nnoreads','nhom_d','nhet_d','nref_d','nnoreads_d', 'match','nomatch','matchpct')]

# Now evoke NICE genotypes
gts_full = vcf[[2]]
gts_full = cut_cm_to_samples(gts_full, used_samples_uniq, 5, include_children)
gts_full_interesting = gts_full[!names(gts_full) %in% c('REF','ID','ALT','QUAL','FILTER','INFO','FORMAT')]
gts_full_interesting <- as.data.frame(lapply(gts_full_interesting, as.character))
cm_final = left_join(cm5, gts_full_interesting, by=c('chrom','start'))
colnames(cm_final)

#optional
sourceding = cmah[,c('chrom','start','end','CALLERSET_LIST')]
sourceding <- as.data.frame(lapply(sourceding, as.character))

cm_final = left_join(sourceding, cm_final, by=c('chrom','start','end'))

outresname = 'review_peter4xx'

# write.table(cm_final, file=paste0('~/PhD/projects/huminvs/mosaicatcher/analysis/results/',outresname,'/arbigent_results/reGTs.csv'), sep='\t',col.names = T,row.names=F, quote = F)
# write.table(cm_final[cm_final$verdict=='MendelFail',], file=paste0('~/PhD/projects/huminvs/mosaicatcher/analysis/results/',outresname,'/arbigent_results/reGTs_MendelFail.csv'), sep='\t',col.names = T,row.names=F, quote = F)
# write.table(cm_final[cm_final$verdict=='NoEvents',], file=paste0('~/PhD/projects/huminvs/mosaicatcher/analysis/results/',outresname,'/arbigent_results/reGTs_NoEvents.csv'), sep='\t',col.names = T,row.names=F, quote = F)
# write.table(cm_final[cm_final$verdict=='AlwaysComplex',], file=paste0('~/PhD/projects/huminvs/mosaicatcher/analysis/results/',outresname,'/arbigent_results/reGTs_AlwaysComplex.csv'), sep='\t',col.names = T,row.names=F, quote = F)
# # 
# write.table(cm_final[!(cm_final$verdict %in% c('MendelFail','NoEvents','AlwaysComplex')),], file=paste0('~/PhD/projects/huminvs/mosaicatcher/analysis/results/',outresname,'/arbigent_results/reGTs_survivors.csv'), sep='\t',col.names = T,row.names=F, quote = F)
# # 
# 
# gts_nice = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/review_peter4xx/arbigent_results/U32_freezemerge.csv'
# gtn = read.table(gts_nice, header=T, sep='\t')
# gtn = cut_cm_to_samples(gtn, used_samples_uniq, 7, include_children)
# gtn <- as.data.frame(lapply(gtn, as.character))
# cm_final2 = left_join(cm5,gtn, by=c('chrom','start','end'))
# cm_final2_cut = cm_final2[,!(names(cm_final2) %in% c('ID.y','len.y','valid_bins.y'))]
# cm_final2_cut = left_join(sourceding, cm_final2_cut, by=c('chrom','start','end'))
# 
# write.table(cm_final2_cut, file=paste0('~/PhD/projects/huminvs/mosaicatcher/analysis/results/',outresname,'/arbigent_results/reGTs_more_readable.csv'), sep='\t',col.names = T,row.names=F, quote = F)
# write.table(cm_final2_cut[!(cm_final2_cut$verdict %in% c('MendelFail','NoEvents','AlwaysComplex')),], file=paste0('~/PhD/projects/huminvs/mosaicatcher/analysis/results/',outresname,'/arbigent_results/reGTs_survivors_more_readable.csv'), sep='\t',col.names = T,row.names=F, quote = F)
#   #                        

######################################
### Here, we will compare phases.  ###
######################################


davidphased = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/bed_factory/revision/david_n35_323/redund_inversions_n35_phased_PAVsync.bed'
minephased = '~/PhD/projects/huminvs/mosaicatcher/analysis/results/review_david323/sv_probabilities/all_rephased.txt'

dp = read.table(davidphased, sep='\t', header=T)
mp = read.table(minephased, sep=' ', header=T)

# Change naming in mine
mp$sample = gsub("GM", "NA", mp$sample)
# Also that HG0002
mp$sample[mp$sample=='HG002'] = 'NA24385'

colnames(mp)[1] = 'seqnames'

both = left_join(mp, dp, by=c('seqnames','start','end','sample'))

# 
# 
# 
# 
# 
# 
# hu = '~/Desktop/oct30/U32_het_mendel.bed'
# hut = read.table(hu, header=1, sep=' ')
# 
# cmah = left_join(cma, hut[,c('chrom','start','end','Major_allele','Major_allele_freq','Heterozygosity','Mendel_cons')])
# cmah = cmah[,c('chrom','start','end','ID','len','valid_bins','mapability','PercOverlap_SD98', 'type','WSSD.med.cn', 'CR','eventsclaimed','match','total','nhom','nhet','nref','nnoreads','ncomplex','CALLERSET_LIST','Major_allele','Major_allele_freq','Heterozygosity','Mendel_cons',used_samples_uniq)]
# write.table(cmah, file='~/Desktop/oct30/bigtable_post.csv', quote = F, col.names = T, row.names=F, sep='\t')
# 
# # #tsne because i'm bored
# # library(Rtsne)
# trn = na.omit(data.matrix(cma[,c('mapability','CR','nhom','nhet','nref','ncomplex', 'eventsclaimed')]))
# trn = na.omit(data.matrix(cma[,c('mapability','CR')]))
# 
# trn = trn/colSums(trn)
# # tsne <- Rtsne(as.matrix(na.omit(trn)), check_duplicates = FALSE, pca = FALSE, perplexity=30, theta=0.5, dims=2)
# # 
# # #library(devtools)
# # #install_github('sinhrks/ggfortify')
# # #library(ggfortify); library(ggplot2)
# # pca_res = prcomp(trn, scale. = TRUE)
# # autoplot(pca_res)
# # 
# # ggplot() + geom_point(aes(x=pca_res$x[,'PC1'], pca_res$x[,'PC2'], color= na.omit(cma[,c('mapability','CR','nhom','nhet','nref','ncomplex','CALLERSET_LIST')]$CALLERSET_LIST)))
# # 
# # library(devtools)
# # install_github('sinhrks/ggfortify')
# # library(ggfortify); library(ggplot2)
# # data(iris)
# # iris.pca <- iris[c(1, 2, 3, 4)] 
# # autoplot(prcomp(iris.pca))
# 
# 
# #2*rowSums2(cm=='1|1')) + rowSums2(cm=='1|0') + rowSums2(cm=='0|1') + 
# #   (2*rowSums2(cm=='1|1_lowconf')) + rowSums2(cm=='1|0_lowconf') + rowSums2(cm=='0|1_lowconf')
# # cm$n = 2*((rowSums2(cm=='1|1')) + rowSums2(cm=='1|0') + rowSums2(cm=='0|1') + 
# #             (rowSums2(cm=='1|1_lowconf')) + rowSums2(cm=='1|0_lowconf') + rowSums2(cm=='0|1_lowconf') +
# #             (rowSums2(cm=='0|0_lowconf')) + (rowSums2(cm=='0|0')))
# 
# #3D plot
# library("plot3D")
# cmah2 = na.omit(cmah)
# x = cmah2$CR
# y = cmah2$valid_bins
# z = cmah2$Heterozygosity
# col = z
# scatter3D(x,y,z,colval=z, col=NULL, add=FALSE)
# 
# ggplot(cmah2) + geom_point(aes(x=CR, y=Heterozygosity, color=Heterozygosity))# +
# scale_x_log10() + scale_y_log10()
# 
# aa = cmah2[cmah2$valid_bins<10,]


ggplot(cm_final) + geom_bar(aes(x=verdict))
