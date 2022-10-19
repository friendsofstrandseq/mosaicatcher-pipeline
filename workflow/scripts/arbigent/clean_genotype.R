# Genotype
library(ggplot2)
library(reshape)
library(dplyr)
library(tibble)
source('clean_genotype_defs.R')
source('vcf_defs.R')
# but first explore

runname = 'review_david323'
debug_file = paste0('~/PhD/projects/huminvs/mosaicatcher/analysis/results/', runname, '/msc.debug')
outdir = paste0('/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/', runname, '/arbigent_results')



load_tab <- function(){
  #full_list = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/peter_all_oct16/sv_probabilities/all.txt'
  #full_list = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/Sseq25/sv_probabilities/all.txt'
  #full_list = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/David_55overhang/sv_probabilities/all.txt'
  #full_list = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/HG00733_oct12/res2/all.txt'
  #full_list = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/U32_freezemerge2/sv_probabilities/all_rephased.txt'
  #verdicts_link = '~/PhD/projects/huminvs/mosaicatcher/david_list_15sep/our_verdicts_oct1.txt'
  #full_list = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/U32_freezemerge2/'
  full_list = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/review_david323/sv_probabilities/all.txt'
  
  
  
  david_list = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/bed_factory/merge_27oct/sv_inv_bu.bed'
  dl = read.table(david_list, header=T, sep='\t')
  dli = dl[c('CHROM','POS','END','ID','MERGE_SAMPLES', 'CALLERSET_LIST', 'PercOverlap_SD98', 'type', 'WSSD.med.cn')]
  
  tab = read.table(full_list, header=T, stringsAsFactors = F)
  #verd = read.table(verdicts_link, header=T, stringsAsFactors = F)
  #verd = verd[verd$chr!='chrY',]
  
  #tab = merge(tab, verd)
  tab = tab %>% mutate(ID = paste0(chrom,'-',start+1,'-INV-',(end-start)+1))
  
  return(tab)
}


tab = load_tab()

cutoff = 3 # Cutoff for 'lowconf'

tab[tab=='nomappability']='./.'
# Complex vs simple: we want complex LLHs to be double the ones
# of simple, and at the same time at least higher by magnitude 5.
bias_factor = tab$confidence_hard_over_second # This is the first criterion: double 
bias_add_factor = 5 # Second criterion: plus 5
mpl = 1./bias_factor # Just cosmetics. Not sure if needed even

# tabp: the simple list
tabp = add_gts_revisited_lowconf(tab, bias_factor, bias_add_factor, cutoff)
tabp3 = add_gts_revisited(tab, bias_factor, bias_add_factor, cutoff)

# tabp2: the list including everything
tabp2 = add_long_gts_revisited(tab, bias_factor, bias_add_factor, cutoff)

# simple
#tab = add_simple_gts(tabp)

### PLOT A BIT ### 
tab2 = tab[tab$sample=='HG00733',]
# ggplot(data=tab) + geom_point(aes(x=(confidence_hard_over_second+1), 
#                                   y=(confidence_nobias_over_hard+
#                                        confidence_hard_over_second+1), 
#                                   color=GTs), size=1) + 
#   scale_x_log10() + 
#   scale_y_log10() + 
#   coord_fixed(ratio = 1,expand = TRUE) +
#   labs(x='Confidence simple inversion [LLH inv / LLH secondhighestinv]', 
#        y='Most confident event [LLH event / secondhighestinv]',
#        title='Prediction likelihoods')

CN = read.table(debug_file, header=1, stringsAsFactors = F)
CNmerge = as.data.frame(lapply(CN[, c("chrom","start","end","valid_bins")], as.character))

CNmerge = as.tbl(CNmerge)
CNmerge <- CNmerge %>%  mutate(chrom = as.character(chrom),
                               start = as.numeric(as.character(start)),
                               end = as.numeric(as.character(end)))#,
#                               CN = as.numeric(as.character(CN)))#,
#mapability = as.numeric(as.character(mapability)))
tab3 <- full_join(tab, CNmerge, by = c("chrom","start","end"))
tab3$valid_bins = as.numeric(as.character(tab3$valid_bins))
tabp3 <- full_join(tabp3, CNmerge, by = c("chrom","start","end"))
tabp3$valid_bins = as.numeric(as.character(tabp3$valid_bins))
tabp <- full_join(tabp, CNmerge, by = c("chrom","start","end"))
tabp$valid_bins = as.numeric(as.character(tabp$valid_bins))
tab4 = tab3[tab3$sample=='HG00733',]


callmatrixpass = cast(tab[tab$verdict=='pass',], chr+start+end+width~sample, value='GT')

callmatrix = cast(unique(tabp3), chrom+start+end+ID+len+valid_bins~sample, value='GT')

cms = callmatrix
samplenames = colnames(cms)[7:length(colnames(cms))]
cms_part = cms[,samplenames]
cms_part[] <- lapply(cms_part, as.character)
#cm_verysimple[,samplenames][cm_verysimple[,samplenames] %in% '0|1',]
unique(as.vector(as.matrix(cms_part)))
# Strongly simplify things

cms_part[cms_part == 'noreads'] <- './.'
cms_part[cms_part == '0000'] <- './.'
cms_part[cms_part == '2200'] <- './.'
cms_part[cms_part == '2110'] <- './.'
cms_part[cms_part == '3100'] <- './.'
cms_part[cms_part == '1101'] <- './1'
cms_part[cms_part == '0100'] <- '1|.'
cms_part[cms_part == '2101'] <- './1'
cms_part[cms_part == '1110'] <- './0'
cms_part[cms_part == '0010'] <- './1'
cms_part[cms_part == '1120'] <- './.'
cms_part[cms_part == '1210'] <- './.'
cms_part[cms_part == '0103'] <- '1|1'
cms_part[cms_part == '1201'] <- './1'
cms_part[cms_part == '0001'] <- './1'
cms_part[cms_part == '1102'] <- './.'
cms_part[cms_part == '1030'] <- '0|0'
cms_part[cms_part == '1200'] <- './.'
cms_part[cms_part == '0301'] <- '1|1'
cms_part[cms_part == '2000'] <- './.'
cms_part[cms_part == '4000'] <- './.'
cms_part[cms_part == '3010'] <- '0|0'
cms_part[cms_part == '0120'] <- '1|0'
cms_part[cms_part == '1100'] <- './.'
cms_part[cms_part == '3001'] <- './1'
cms_part[cms_part == '0030'] <- './.'
cms_part[cms_part == '2020'] <- '0|0'
cms_part[cms_part == '0202'] <- '1|1'
cms_part[cms_part == '1300'] <- './.'
cms_part[cms_part == '1000'] <- './.'
cms_part[cms_part == '2100'] <- './.'

cms_full = cbind(cms[,1:6], cms_part)


callmatrix = cast(unique(tabp), chrom+start+end+ID+len+valid_bins~sample, value='GT')

#callmatrix = cast(tabp3, chrom+start+end+ID+len+valid_bins~sample, value='GT')

#write.table(callmatrix, file='~/Desktop/oct21/U32_freezemerge.csv', quote = F, col.names = T, row.names=F, sep='\t')
callmatrix_detail = cast(tabp, chrom+start~sample, value='GT')


callmatrix_detail = cast(tabp2, chrom+start+end+ID+len+orig_label~sample, value='GTL')
#write.table(callmatrix_detail, file='~/Desktop/oct21/U32_freezemerge_detail.csv', quote = F, col.names = T, row.names=F, sep='\t')
callmatrix_detail = callmatrix_detail[ , colSums(is.na(callmatrix_detail)) == 0]

# Sidequest: find hom invs
callmatrix_hom_lab = callmatrix
callmatrix_hom_lab$nhom = rowSums(callmatrix_hom_lab=='1|1')
cm_hom = callmatrix_hom_lab[callmatrix_hom_lab$nhom > (dim(callmatrix_hom_lab)[2] - 5)*0.8,]
hom_ids = cm_hom$ID
cm_detail_hom = callmatrix_detail[callmatrix_detail$ID %in% hom_ids,]


vcf = vcfify_callmatrix_detail(callmatrix_detail)
vcf_miso = vcfify_callmatrix_detail(cm_detail_hom)
vcf_limix = vcfify_callmatrix_simple_for_limix(cms_full)
save=F
if (save==T){
  
  # Prep directory
  dir.create(outdir)
  
  # Paths, paths, paths. 
  callmatrix_file = 'U32_freezemerge.csv'
  callmatrix_file_detail = 'U32_freezemerge_detail.csv'
  vcffile_all = 'res_all.vcf'
  vcffile_miso = 'res_miso.vcf'
  vcffile_limix = 'res_verysimple.vcf'
  
  # Save simple callmatrix
  write.table(callmatrix, file=file.path(outdir,callmatrix_file), quote = F, col.names = T, row.names=F, sep='\t')
  # Same detailed callmatrix
  write.table(callmatrix_detail, file=file.path(outdir, callmatrix_file_detail), quote = F, col.names = T, row.names=F, sep='\t')
  # Save vcf_all
  outvcffile = file.path(outdir, vcffile_all)
  writeLines(vcf[[1]], file(outvcffile))
  write.table(vcf[[2]], file=outvcffile, quote=F, col.names=F, row.names=F, sep='\t', append = T)
  system(paste0('bcftools sort -O v -o ', outvcffile, ' ', outvcffile))
    
  # Save vcf_miso
  outvcffile_miso = file.path(outdir, vcffile_miso)
  writeLines(vcf_miso[[1]], file(outvcffile_miso))
  write.table(vcf_miso[[2]], file=outvcffile_miso, quote=F, col.names=F, row.names=F, sep='\t', append = T)
  system(paste0('bcftools sort -O v -o ', outvcffile_miso, ' ', outvcffile_miso))
  
  # Save vcf_limix
  outvcffile_limix = file.path(outdir, vcffile_limix)
  writeLines(vcf_limix[[1]], file(outvcffile_limix))
  write.table(vcf_limix[[2]], file=outvcffile_limix, quote=F, col.names=F, row.names=F, sep='\t', append = T)
  system(paste0('bcftools sort -O v -o ', outvcffile_limix, ' ', outvcffile_limix))
  
  # Save only complex stuff Complex exploration
  tabcomp = tab[tab$GTs == 'complex',]
  indiv_invs <- unique( tabcomp[ , 1:3 ] )
  library(dplyr)
  aa = left_join(indiv_invs,callmatrix_detail)
  bb = aa[
    with(aa, order(chrom, start)),
    ]
  
  #write.table(bb, file='~/Desktop/oct21/complex_calls.csv', quote = F, col.names = T, row.names=F, sep='\t')
  }



# Sidequest: find hom invs
callmatrix2 = cast(tabp3, chrom+start+end+ID+len~sample, value='GT')

callmatrix_hom_lab = callmatrix2
callmatrix_hom_lab$nhom = rowSums(callmatrix_hom_lab == '1|1')
cm_hom = callmatrix_hom_lab[callmatrix_hom_lab$nhom > (dim(callmatrix_hom_lab)[2] - 5)*0.9,]
hom_ids = cm_hom$ID
cm_detail_hom = callmatrix_detail[callmatrix_detail$ID %in% hom_ids,]

