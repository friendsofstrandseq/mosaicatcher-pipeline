# Whoeps, 30th Nov 2020
# Flip vcf labels!

library(stringr)
library(ggplot2)
library(optparse)
library(data.table)
library(dplyr)

make_nophasing_list <- function(phases_proc_rename_f, blacklist_f, blacklist_lower, blacklist_higher){
  #print(phases_proc_rename_f)
  samples_no_phasing = phases_proc_rename_f
  samples_no_phasing$diagnose = 'Sync_successful'
  print('ho') 
  for (sample in samples_no_phasing$sample){
      if (all(samples_no_phasing[samples_no_phasing$sample==sample,]$pct_match == 100)){
	    samples_no_phasing[samples_no_phasing$sample==sample,]$diagnose = "PAV_vcf_missing"
     }
	  if (all(samples_no_phasing[samples_no_phasing$sample==sample,]$pct_match == 0)){
	    samples_no_phasing[samples_no_phasing$sample==sample,]$diagnose = "No_phasing_required"
	 }
   }
  #samples_no_phasing[(samples_no_phasing$pct_match == 100),]$diagnose = 'PAV_vcf_missing'
  #samples_no_phasing[(samples_no_phasing$pct_match == 0),]$diagnose = 'No_phasing_required'
  
  # Blacklist samples automatically. 
  print(blacklist_lower)
  print(blacklist_higher)
  samples_no_phasing[((samples_no_phasing$pct_match > blacklist_lower) & (samples_no_phasing$pct_match < blacklist_higher)),]$diagnose = 'Blacklisted'
  print('bye')
  samples_no_phasing_report = samples_no_phasing[,c('sample','sseqsample','chr','pct_match', 'diagnose', 'flip')]
return(samples_no_phasing_report)
}

beautify_phases <- function(phases_f){
  colnames(phases_f) = c('sample_chr_path', 'V2', 'V3',
                         'sample', 'V5', 'V6', 'pct_match',
                         'V8','V9','V10')
  
  phases_f$chr = str_split_fixed(phases_f$sample_chr_path, '/', 10)[,3]
  phases_f$sample = str_split_fixed(phases_f$sample_chr_path, '/', 10)[,2]
  #phases_f$sample = str_split_fixed(phases_f$sample_long, '_', 10)[,2]
  phases_f$flip = phases_f$pct_match < 50
  
  # If there is nothing, please don't invert we have to examine first. 
  # Probably it's a X or Y chromsome
  phases_f[phases_f$pct_match == 0,]$flip = FALSE
  
  return(phases_f[,c('sample','chr','pct_match','flip')])
}

flip_gt <- function(gt_factors){
  
  gt_list = as.character(gt_factors)
  
  # Make sure there are only simple and complex preds (n=3 and n=4)
  stopifnot(length(gt_list[nchar(as.character(gt_list)) == 3]) +
              length(gt_list[nchar(as.character(gt_list)) == 4]) ==
              length(gt_list)
  )
  
  # Flip simple preds
  gt_list[nchar(as.character(gt_list)) == 3] = 
    paste0(substr(gt_list[nchar(as.character(gt_list)) == 3], 3, 3),
           '|',
           substr(gt_list[nchar(as.character(gt_list)) == 3], 1, 1))
  
  # Flip complex preds
  gt_list[nchar(as.character(gt_list)) == 4] = 
    paste0(substr(gt_list[nchar(as.character(gt_list)) == 4], 3, 4),
           substr(gt_list[nchar(as.character(gt_list)) == 4], 1, 2))
  
  return(gt_list)
}

unphase <- function(gt_factors){
  
  gt_list = as.character(gt_factors)
  
  # Make sure there are only simple and complex preds (n=3 and n=4)
  stopifnot(length(gt_list[nchar(as.character(gt_list)) == 3]) +
              length(gt_list[nchar(as.character(gt_list)) == 4]) ==
              length(gt_list)
  )
  
  # Unphase
  gt_list[nchar(as.character(gt_list)) == 3] = 
    paste0(substr(gt_list[nchar(as.character(gt_list)) == 3], 3, 3),
           '/',
           substr(gt_list[nchar(as.character(gt_list)) == 3], 1, 1))
  
  # TODO
  # Complex preds we can't really do anything about them I'm afraid. 
  
  return(gt_list)
}


rename_samples <- function(phases_f, rename_f){
  
  rnames = read.table(rename_f, header=T, sep=',', stringsAsFactors=F)
  
  phases_f$sseqsample = 'unk'
  for(i in 1:nrow(phases_f)) {
    phases_f[i,]$sseqsample =  rnames[rnames$sseq==phases_f[i,]$sample,]$sseqshort
  }
  
  
  return(phases_f)
}

#INPUT INSTRUCTIONS
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="vcf file to be re-phased", metavar="character"),
  make_option(c("-b", "--bed"), type="character", default=NULL,
              help="Table with re-phasing info, created by invphasing snakemake", metavar="character"),
  make_option(c("-s", "--sampletransitions"), type="character", default=NULL,
              help="Table with naming translation table: long sseq to short sseq", metavar="character"),
  make_option(c("-x", "--blacklist"), type="character", default=NULL,
              help="Tab-separated blacklist [sample\tchr] of chrs not to phase", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="./outputcorr/",
              help="Outputdir for phased all.txt and other qcs", metavar="character")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

regenotyper_res = opt$file
phases_file = opt$bed
renaming_file = opt$sampletransitions
blacklist_file = opt$blacklist
outdir = opt$outdir
blacklist_lower = 20
blacklist_higher = 80
# regenotyper_res = '~/scratch/invphasing/input_regenotyper/all.txt'
# phases_file = '/home/hoeps/scratch/invphasing/res_backup2/similaritymerge.txt'
# #phases_file ='~/scratch/invphasing/res/similaritymerge.txt'
# renaming_file = '~/scratch/invphasing/input_naming/samplenames_sseq_to_sseq_short.txt'
# blacklist_file = '~/scratch/invphasing/blacklist/sample_chr_blacklist.txt'
# outdir = '~/Desktop/lab/lab2'

#regenotyper_res = "~/s/g/korbel/hoeps/projects/huminvs/mosai_results/results_hg38_run4filter/regenotyper_allsamples_bulk/all_sv_calls_unphased.txt"
#phases_file = '~/s/g/korbel2/StrandSeq/Test_WH/pipeline_7may/pipeline/utils/regenotyper/phasing_res/similaritymerge.txt'
#blacklist_file = '~/s/g/korbel2/StrandSeq/Test_WH/pipeline_7may/pipeline/utils/regenotyper/blacklist/sample_chr_blacklist.txt'
#renaming_file = '~/s/g/korbel2/StrandSeq/Test_WH/pipeline_7may/pipeline/utils/regenotyper/input_naming/samplenames_sseq_to_sseq_short.txt'
#outdir = '~/Desktop/labphase'


all = read.table(regenotyper_res, sep=' ', header=T, stringsAsFactors = F)
allbackup = all 

# Prepare phases: beautify and rename
phases = read.table(phases_file, header=0, stringsAsFactors = F)
phases_proc = beautify_phases(phases)
phases_proc_rename = rename_samples(phases_proc, renaming_file)


# keep track of what we have seen:
all_samples = unique(all$sample)
samplecheck_df = data.frame(sample = all_samples, seen = 'F', stringsAsFactors = F)
samplecheck_df[samplecheck_df$sample %in% unique(phases_proc_rename$sseqsample),]$seen = 'T'

# blacklist
blacklist = read.table(blacklist_file, header=T, sep='\t')



# tie together sample and chr for easier df subsettig. 
phases_proc_rename$sample_chr = paste0(phases_proc_rename$sseqsample, phases_proc_rename$chr)
blacklist$sample_chr = paste0(blacklist$sseqsample, blacklist$chr)

# Prevent blacklisted from flipping
if ((dim(blacklist)[1]>0) && (any(phases_proc_rename$sample_chr %in% blacklist$sample_chr))){
  phases_proc_rename[phases_proc_rename$sample_chr %in% blacklist$sample_chr,]$flip = 'Blacklisted'
}

print('got so far')
# keep track of non-valid PAV samples.
samples_no_phasing_report = make_nophasing_list(phases_proc_rename, blacklist, blacklist_lower, blacklist_higher)

print(dim(phases_proc_rename[samples_no_phasing_report$diagnose == "Blacklisted",])[1])
if (dim(phases_proc_rename[samples_no_phasing_report$diagnose == "Blacklisted",])[1] > 0){
	phases_proc_rename[samples_no_phasing_report$diagnose == "Blacklisted",]$flip = FALSE
	phases_proc[samples_no_phasing_report$diagnose == "Blacklisted",]$flip = FALSE
    samples_no_phasing_report[samples_no_phasing_report$diagnose == "Blacklisted",]$flip = FALSE
}


# all.txt can contain 'nomappability' entries. If that is the case, replace them with ./.
if ('nomappability' %in% all$pred_hard){
  all[all=='nomappability'] = './.'
}


# Now the fuzzy part
# loop over samples and chrs of phases_proc. Unseen samples are never considered here (!)

# [W] Update 11th March: Actually, the Blacklisted samples we also want to keep with their SSeq phasing. 
# So we're not 'unphasing' them (commented out as you can see), and report later w.r.t what everything
# is phased.  
for (sample in unique(phases_proc_rename$sseqsample)){
  
  for (chr in unique(phases_proc_rename$chr)){
    # Do GTs on this chr need to be flipped?
    # if needs to flip (flip = True)
    if (phases_proc_rename[(phases_proc_rename$sseqsample==sample) & (phases_proc_rename$chr==chr),]$flip == T){
      #print(paste0('Flipping ', sample, ' ', chr))
      all[all$sample==sample & all$chrom==chr,]$pred_hard = flip_gt(all[all$sample==sample & all$chrom==chr,]$pred_hard)
      all[all$sample==sample & all$chrom==chr,]$pred_soft = flip_gt(all[all$sample==sample & all$chrom==chr,]$pred_soft)
      all[all$sample==sample & all$chrom==chr,]$pred_nobias = flip_gt(all[all$sample==sample & all$chrom==chr,]$pred_nobias)
      all[all$sample==sample & all$chrom==chr,]$second_hard = flip_gt(all[all$sample==sample & all$chrom==chr,]$second_hard)
    } #else if (phases_proc_rename[(phases_proc_rename$sseqsample==sample) & (phases_proc_rename$chr==chr),]$flip == 'Blacklisted'){
    #print(paste0('Blacklisted: ', sample, ' ', chr))
    #all[all$sample==sample & all$chrom==chr,]$pred_hard = unphase(all[all$sample==sample & all$chrom==chr,]$pred_hard)
    #all[all$sample==sample & all$chrom==chr,]$pred_soft = unphase(all[all$sample==sample & all$chrom==chr,]$pred_soft)
    #all[all$sample==sample & all$chrom==chr,]$pred_nobias = unphase(all[all$sample==sample & all$chrom==chr,]$pred_nobias)
    #all[all$sample==sample & all$chrom==chr,]$second_hard = unphase(all[all$sample==sample & all$chrom==chr,]$second_hard)
  #else {
    #print(paste0('No Flipping ', sample, ' ', chr))
  #}
}
}

# Additionally, unphase all calls in samples which do not appear in the phases file. Ideally, this should be none. 
unprocessed_samples = samplecheck_df[samplecheck_df$seen==F,]$sample
all[(all$sample %in% unprocessed_samples),]$pred_hard = unphase(all[(all$sample %in% unprocessed_samples),]$pred_hard)
all[(all$sample %in% unprocessed_samples),]$pred_soft = unphase(all[(all$sample %in% unprocessed_samples),]$pred_soft)
all[(all$sample %in% unprocessed_samples),]$pred_nobias = unphase(all[(all$sample %in% unprocessed_samples),]$pred_nobias)
all[(all$sample %in% unprocessed_samples),]$second_hard = unphase(all[(all$sample %in% unprocessed_samples),]$second_hard)


# make/write LOG
samples_with_pav = samples_no_phasing_report[samples_no_phasing_report$diagnose == 'Sync_successful',]

l1 = '## GT flipping log ##'
l2 = paste0('Samples in regenotyper input: ', length(unique(allbackup$sample)))
l3 = paste0('Samples in phasing table: ', length(unique(phases_proc$sample)))
l3.1 = paste0('Samples without PAV: ', length(unique(samples_no_phasing_report[samples_no_phasing_report$diagnose=='PAV_vcf_missing',]$sample)))
l3.2 = paste0('...which are: ')
l3.21 = as.character(unique(samples_no_phasing_report[samples_no_phasing_report$diagnose=='PAV_vcf_missing',]$sample))
l3.3 = paste0('Chrs Blacklisted: ', dim(unique(samples_no_phasing_report[samples_no_phasing_report$diagnose=='Blacklisted',c('sample', 'chr')]))[1])
l3.4 = paste0('...which are: ')
l3.41 = as.character(unique(samples_no_phasing_report[samples_no_phasing_report$diagnose=='Blacklisted',c('sample', 'chr')]))
l4 = paste0('Haplotype labels changed in ', sum(phases_proc$flip) ,' out of ', length(phases_proc$flip), ' chromosomes considered (', round(sum(phases_proc$flip)/length(phases_proc$flip)*100,0),'%)')
l5 = paste0('If only considering samples with PAV:')
l6 = paste0('Haplotype labels changed in ', sum(as.logical(samples_with_pav$flip)) ,' out of ', length(samples_with_pav$flip), ' chromosomes considered (', round(sum(as.logical(samples_with_pav$flip))/length(samples_with_pav$flip)*100,0),'%)')
log_all = c(l1,l2,l3,l3.1,l3.2,l3.21,l3.3,l3.4,l3.41,l4,l5,l6,collapse='\n', sep='\n')




outlogfile = file.path(outdir, 'phasing.log')
writeLines(log_all, file(outlogfile))

# save nophasing report
outnophasereportfile = file.path(outdir, 'rephase_chr_log.txt')
write.table(samples_no_phasing_report, outnophasereportfile, row.names=F, col.names=T, quote=F, sep='\t')

# save samplecheck_df
outsamplecheckfile = file.path(outdir, 'samplecheck')
write.table(samplecheck_df, outsamplecheckfile, row.names=F, col.names=T, quote=F, sep='\t')


cols = rep(c('#3a32c9','#87c1ff'),100)
# make/save a QC plot
p = ggplot(samples_no_phasing_report[samples_no_phasing_report$diagnose != 'PAV_vcf_missing',]) + 
  geom_point(aes(x=paste(sample,chr), y=pct_match, color=sample), size = 0.5) + 
  #theme(axis.text.x = element_text(angle = 90)) + 
  labs(title = 'Percent identity phased snps: Strandseq vs PAV phasing', x= 'Sample and Chr', y= 'Identity') +
  theme(#axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        legend.title = element_blank(),
        legend.position = "none"
        ) + 
  scale_color_manual(values=cols) + 
  geom_vline(aes(xintercept=0)) + 
  geom_hline(aes(yintercept=0)) + 
  geom_hline(aes(yintercept=100)) + 
  geom_abline(aes(slope=0, intercept=20), size=0.1) + 
  geom_abline(aes(slope=0, intercept=80), size=0.1)
#p
outplotfile = file.path(outdir, 'phasing_qc.pdf')
ggsave(filename=outplotfile, plot=p, width=15, height=6.66, units='cm', device='pdf')

# save phase_f
outphasefile = file.path(outdir, 'phasing_clean.txt')
write.table(phases_proc_rename, file=outphasefile, row.names=F, col.names = T, quote = F, sep='\t')

# save all_phased
outphasefile = file.path(outdir, 'all_rephased.txt')
write.table(all, file=outphasefile, row.names=F, col.names = T, quote = F, sep=' ')



