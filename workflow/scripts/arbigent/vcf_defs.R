# Whoeps, 23 Oct 2020
# Function(s) needed to create a vcf file that will be saved as a final result of arbigent



#' Take a callmatrix and make a vcf out of it. 
#' @param callmatrix
#' @return list: [all header lines, all data lines]
#' @author Wolfram Hoeps
#' @export
vcfify_callmatrix_detail <- function(cd){
  # Header lines
  l1 = '##fileformat=VCFv4.2'
  l2 = paste0('##fileDate=',Sys.Date())
  l3 = '##ALT=<ID=INV,Description="Inversion">'
  l4 = '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">'
  l5 = '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">' 
  l6 = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Main Genotype">'
  l62 = '##FORMAT=<ID=CF,Number=1,Type=String,Description="High Confidence (T/F) in Main Genotype?">'
  l7 = '##FORMAT=<ID=AT,Number=1,Type=String,Description="Alternative Genotype in case Main is complex">'
  l8 = '##FORMAT=<ID=GL,Number=G,Type=String,Description="Log10-scaled genotype likelihoods for Main Genotype">'
  l9 = '##FORMAT=<ID=AL,Number=G,Type=String,Description="Log10-scaled genotype likelihoods for Alternative Genotype">'
  l10 = '##reference=/g/solexa/bin/genomesNew/GRCh38Decoy/GRCh38Decoy.fa'
  l11 = '##contig=<ID=chr1,length=248956422>'
  l12 = '##contig=<ID=chr1,length=248956422>'
  l13 = '##contig=<ID=chr2,length=242193529>'
  l14 = '##contig=<ID=chr3,length=198295559>'
  l15 = '##contig=<ID=chr4,length=190214555>'
  l16 = '##contig=<ID=chr5,length=181538259>'
  l17 = '##contig=<ID=chr6,length=170805979>'
  l18 = '##contig=<ID=chr7,length=159345973>'
  l19 = '##contig=<ID=chr8,length=145138636>'
  l20 = '##contig=<ID=chr9,length=138394717>'
  l21 = '##contig=<ID=chr10,length=133797422>'
  l22 = '##contig=<ID=chr11,length=135086622>'
  l23 = '##contig=<ID=chr12,length=133275309>'
  l24 = '##contig=<ID=chr13,length=114364328>'
  l25 = '##contig=<ID=chr14,length=107043718>'
  l26 = '##contig=<ID=chr15,length=101991189>'
  l27 = '##contig=<ID=chr16,length=90338345>'
  l28 = '##contig=<ID=chr17,length=83257441>'
  l29 = '##contig=<ID=chr18,length=80373285>'
  l30 = '##contig=<ID=chr19,length=58617616>'
  l31 = '##contig=<ID=chr20,length=64444167>'
  l32 = '##contig=<ID=chr21,length=46709983>'
  l33 = '##contig=<ID=chr22,length=50818468>'
  l34 = '##contig=<ID=chrX,length=156040895>'
  l35 = '##contig=<ID=chrY,length=57227415>'
  samplenames = colnames(cd)[7:length(colnames(cd))]
  sn=paste(samplenames, collapse='\t')
  header_presamples=paste('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sep="\t")
  header = paste(header_presamples, sn, sep='\t')
  
  header_all = paste(l1,l2,l3,l4,l5,l6,l62,l7,l8,l9,l10,
                     l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,
                     l21,l22,l23,l24,l25,l26,l27,l28,l29,l30,
                     l31,l32,l33,l34,l35,header,collapse='\n', sep='\n')
  
  cprint = cd
  
  cprint = cprint %>% mutate(REF = '.')
  cprint = cprint %>% mutate(ALT = '<INV>')
  cprint = cprint %>% mutate(QUAL = '.')
  cprint = cprint %>% mutate(FILTER = 'PASS')
  cprint = cprint %>% mutate(INFO = paste0('END=', end))
  cprint = cprint %>% mutate(FORMAT = 'GT:CF:AT:GL:AL')

  
  cprint <- data.frame(lapply(cprint, function(x) {
                      gsub("noreads", "./.:noreads:./.:0:0", x)
                  }))
  
  cprint2 = cprint[,c('chrom','start','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', samplenames)]
  return(list(header_all, cprint2))
}

#' Take a callmatrix and make a vcf out of it. 
#' All very simple for limix. GT only.
#' @param callmatrix
#' @return list: [all header lines, all data lines]
#' @author Wolfram Hoeps
#' @export
vcfify_callmatrix_simple_for_limix <- function(cd){
  
  # Header lines
  l1 = '##fileformat=VCFv4.2'
  l2 = paste0('##fileDate=',Sys.Date())
  l3 = '##ALT=<ID=INV,Description="Inversion">'
  l3.2 = '##INFO=<ID=ID,Number=A,Type=String,Description="Variant IDs per ALT allele.">'
  l4 = '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">'
  l5 = '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">' 
  l6 = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Main Genotype">'
  l7 = '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">'
  l8 = '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">'
  l10 = '##reference=/g/solexa/bin/genomesNew/GRCh38Decoy/GRCh38Decoy.fa'
  l11 = '##contig=<ID=chr1,length=248956422>'
  l12 = '##contig=<ID=chr1,length=248956422>'
  l13 = '##contig=<ID=chr2,length=242193529>'
  l14 = '##contig=<ID=chr3,length=198295559>'
  l15 = '##contig=<ID=chr4,length=190214555>'
  l16 = '##contig=<ID=chr5,length=181538259>'
  l17 = '##contig=<ID=chr6,length=170805979>'
  l18 = '##contig=<ID=chr7,length=159345973>'
  l19 = '##contig=<ID=chr8,length=145138636>'
  l20 = '##contig=<ID=chr9,length=138394717>'
  l21 = '##contig=<ID=chr10,length=133797422>'
  l22 = '##contig=<ID=chr11,length=135086622>'
  l23 = '##contig=<ID=chr12,length=133275309>'
  l24 = '##contig=<ID=chr13,length=114364328>'
  l25 = '##contig=<ID=chr14,length=107043718>'
  l26 = '##contig=<ID=chr15,length=101991189>'
  l27 = '##contig=<ID=chr16,length=90338345>'
  l28 = '##contig=<ID=chr17,length=83257441>'
  l29 = '##contig=<ID=chr18,length=80373285>'
  l30 = '##contig=<ID=chr19,length=58617616>'
  l31 = '##contig=<ID=chr20,length=64444167>'
  l32 = '##contig=<ID=chr21,length=46709983>'
  l33 = '##contig=<ID=chr22,length=50818468>'
  l34 = '##contig=<ID=chrX,length=156040895>'
  l35 = '##contig=<ID=chrY,length=57227415>'
  
  cd[,'NA'] = NULL
  samplenames = colnames(cd)[7:length(colnames(cd))]
  sn=paste(samplenames, collapse='\t')
  header_presamples=paste('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sep="\t")
  header = paste(header_presamples, sn, sep='\t')
  
  header_all = paste(l1,l2,l3,l3.2,l4,l5,l6,l7,l8,l10,
                     l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,
                     l21,l22,l23,l24,l25,l26,l27,l28,l29,l30,
                     l31,l32,l33,l34,l35,header,collapse='\n', sep='\n')
  
  cprint = cd
  
  cprint = cprint %>% mutate(REF = '.')
  cprint = cprint %>% mutate(ALT = '<INV>')
  cprint = cprint %>% mutate(QUAL = '.')
  cprint = cprint %>% mutate(FILTER = 'PASS')
  cprint = cprint %>% mutate(INFO = paste0('ID=',ID, ';END=', end,';AC=1;AN=2'))
  cprint = cprint %>% mutate(FORMAT = 'GT')
  
  
  #cprint <- data.frame(lapply(cprint, function(x) {
  #  gsub("noreads", "./.:noreads:./.:0:0", x)
  #}))
  print(head(cd))
  print(samplenames)
  
  cprint2 = cprint[,c('chrom','start','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', samplenames)]
  return(list(header_all, cprint2))
}
  
