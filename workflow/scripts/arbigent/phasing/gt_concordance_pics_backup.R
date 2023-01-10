# Whoeps 5th Dec 2020
# Horizontal GT concordance plot for both.vcf

library(ggplot2)
library(optparse)

# Read input
#INPUT INSTRUCTIONS
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="vcf file to be checked", metavar="character"),
  make_option(c("-s", "--sname"), type="character", default=NULL,
              help="Samplename: for output filename and title", metavar="character"),
  make_option(c("-c", "--chrname"), type="character", default=NULL,
              help="Chrname: for output filename and title", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="outfile", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

vcf_file = opt$file
samplename = opt$sname
chrname = opt$chrname
outfile = opt$outfile

# Testphase
# both_f = '/home/hoeps/scratch/invphasing/isec_output/HVN3GAFXY_HG02587x02_19s004140-1-1/chr14/isec_output/both.vcf'
# samplename = 'HG02587'
# chr = 'chr1'
#vcf_file = '/home/hoeps/scratch/invphasing/isec_output/NA19239/chrX/isec_output/both.vcf'

# Sometimes this vcf is empty. The tryCatch takes care of that (although it is not an ideal solution. But if it works it works)
tryCatch(
  # If vcf has entries and all cool ...
  expr = {
    both = read.table(vcf_file, stringsAsFactors = F)
    both = both[,c('V1','V2','V10','V11')]
    colnames(both) = c('chr','pos','GTssfull','GTpavfull')
    
    both$GTss = substr(both$GTssfull,1,3)
    both$GTpav = substr(both$GTpavfull,1,3)
    
    both = both[,c('chr','pos','GTss','GTpav')]
    both$pos = both$pos / 1000000
    both$sim = both$GTss == both$GTpav
    both$simconnect = (both$sim-0.5)*2
    
    
    p = ggplot(both) + 
      geom_point(aes(x=pos, y=simconnect, color=sim)) +   
      geom_segment(aes(x=pos, 
                       xend=pos, 
                       y=0, 
                       yend=simconnect,
                       color=sim,
      ),
      size=0.1) +
      labs(title=paste0('HET SNPs: Phase concordance Sseq vs PAV\n', samplename, ' ', chrname), x='Position [Mb]', y='SNP: Phase concordance')
    
    # I want dynamic plot width depending on how long the chr is. 250[Mb] is roughly chr1, so this is our reference.
    # It's all not super exact so yeah, doesnt matter. Also the space that the legend takes isn't taken into account.
    width = 30#(max(both$pos)/250)*50
    ggsave(filename=paste0(outfile,'_smaller.png'), plot=p, width=width, height=9, units='cm', device='png')
    ggsave(filename=outfile, plot=p, width=width, height=9, units='cm', device='pdf')
    },
  # Otherwise, just save some empty plots to satisfy snakemake.
  error = function(e){
    p=ggplot()
    width = 30
    ggsave(filename=paste0(outfile,'_smaller.png'), plot=p, width=width, height=9, units='cm', device='png')
    ggsave(filename=outfile, plot=p, width=width, height=9, units='cm', device='pdf')
    }
)


       