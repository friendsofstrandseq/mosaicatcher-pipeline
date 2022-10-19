# Genotypes QC

# plot 1

library(ggplot2)
library(dplyr)
#theme_set(theme_classic())

#tab2 = tab
tab=tabp
tab=tab[!(tab$pred_hard=='./.'),]
tab$simple = 'ccomplex'
tab[tab$GT %in% c('0|0', '1|0', '0|1', '1|1'),]$simple = 'asimple'
tab[tab$GT %in% c('noreads'),]$simple = 'zeroreads'
tab[tab$GT %in% c('noreads'),]$simple = 'zeroreads'
tab[endsWith(tab$GT, '_lowconf'),]$simple = 'complex_lowconf'
tab[tab$GT %in% c('0|0_lowconf', '1|0_lowconf', '0|1_lowconf', '1|1_lowconf'),]$simple = 'bsimple_lowconf'
tab$verdict = tab$simple
# Histogram on a Continuous (Numeric) Variable
#tab$GT = tab$GTs
g <- ggplot(tab, aes(verdict))# + scale_fill_brewer(palette = "Spectral")

g + geom_histogram(aes(fill=GT), 
                   col="black", 
                   size=.1,
                   stat='count') +   # change number of bins
  labs(title="Genotype predictions for 323 inversions in many samples") 


#### mendel ####

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
ce = make_child_expect_vector()

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
  
  return(F)
  
}


#c2 = callmatrix
callmatrix = c2
callmatrix[callmatrix=='./.']='0|0'
#callmatrix = c2[c2$verdict=='pass',]
callmatrix$mendel1 = 'UNK'
callmatrix$mendel2 = 'UNK'
callmatrix$mendel3 = 'UNK'
callmatrix$mendel4 = 'UNK'

for (row in 1:nrow(callmatrix)){
  #callmatrix[row,]$mendel = 
    callmatrix[row,]$mendel1 = as.logical(test_mendel(ce, callmatrix[row,]$NA19238, callmatrix[row,]$NA19239,callmatrix[row,]$NA19240 ))
    callmatrix[row,]$mendel2 = as.logical(test_mendel(ce, callmatrix[row,]$HG00512, callmatrix[row,]$HG00513,callmatrix[row,]$HG00514 ))
    callmatrix[row,]$mendel3 = as.logical(test_mendel(ce, callmatrix[row,]$HG00731, callmatrix[row,]$HG00732,callmatrix[row,]$HG00733 ))
    callmatrix[row,]$mendel4 = as.logical(test_mendel(ce, callmatrix[row,]$GM19650, callmatrix[row,]$HG00864,callmatrix[row,]$HG03371 ))
    
}
ctest = callmatrix
ct2 = ctest[,c('chr','verdict','HG00731','HG00732','HG00733','HG00512','HG00513','HG00514','NA19238','NA19239','NA19240','mendel1','mendel2','mendel3')]
ct3 = ctest[,c('HG00731','HG00732','HG00733','HG00512','HG00513','HG00514','NA19238','NA19239','NA19240')]

acceptable =  c('0|0', '1|0', '0|1', '1|1',  '0|0_lowconf', '1|0_lowconf', '0|1_lowconf', '1|1_lowconf')
ct4 = ct3[ct3 %in% acceptable,]
#ct3 = ct2[]

#callmatrix[callmatrix$mendel]
#callmatrix$mendelall = callmatrix$mendel1 && callmatrix$mendel2 && callmatrix$mendel3

# 'FALSE' %in% callmatrix[,c('mendel1','mendel2','mendel3')]
# 
mendelall = callmatrix %>% mutate(mm = as.logical(mendel1) * as.logical(mendel2) * as.logical(mendel3))
#mendelall = callmatrix %>% mutate(mm = as.numeric(mendel1) + 1)

# test_mendel(ce, p1, p2, c)
# # First trio

# Histogram on a Continuous (Numeric) Variable
g <- ggplot(mendelall, aes()) + scale_fill_brewer(palette = "Spectral")

g + geom_histogram(aes(fill=as.logical(mm), x=mm), 
                   col="black", 
                   size=.1,
                   stat='count') +   # change number of bins
  labs(title="Genotype predictions for 255 inversions in 31 samples") 

ctest = callmatrix[,c('HG00512','HG00513','HG00514','HG00731','HG00732','HG00733','NA19238','NA19239','NA19240','mendel1','mendel2','mendel3')]
#### compare to david's calls

david_list = '~/PhD/projects/huminvs/mosaicatcher/bed_factory/revision/david_n35_323/nonred_inversions_n35_genotypes.csv'
dgt = read.table(david_list, header=T, sep=',')

#sth else
dq = '~/Desktop/desktop_5th_oct/davidquick.txt'
d = read.table(dq, header=T)
aa = as.data.frame.matrix(table(d$notes, d$phase.pass))
aa$pct = aa[['TRUE']] / (aa[['FALSE']] + aa[['TRUE']])
aa$names = row.names(aa)
aa = aa[order(aa$pct),]
colnames(aa) = c('fail','pass','pct','names')
aa$good = c('bad','nomap','bad','bad','good','bad','good','good','good','bad','bad','good')
dodgewidth <- position_dodge(width=0.9)
ggplot(data=aa, aes(x=names, y=pct), ylim=c(0,1)) + geom_bar(stat='identity', aes(fill=good)) +
  geom_text(aes(label = fail+pass, x = names, y = pct), position = position_dodge(width = 0.8), vjust = -0.6)
  library(pheatmap)
pheatmap(aa, cluster_rows = F)
aa$pct = aa$TRUE

aa3 = aa[,c('names','pct')]
# Histogram on a Continuous (Numeric) Variable
aa2 = melt(aa)
g <- ggplot(aa3, aes(names)) + scale_fill_brewer(palette = "Spectral")

g + geom_histogram(aes(fill=pct), 
                   col="black", 
                   size=.1,
                   stat='count') +   # change number of bins
  labs(title="Genotype predictions for 255 inversions in 31 samples") 

library(ggplot2)
g <- ggplot(d, aes(log10(width)))
g + geom_density(aes(fill=factor(phase.pass)), alpha=0.8) + 
  labs(title="phasing pass vs width", 
       subtitle="City Mileage Grouped by Number of cylinders",
       caption="Source: mpg",
       x="Inv width",
       fill="phase pass")

g <- ggplot(d, aes(log10(width)))
g + geom_density(aes(fill=factor(phase.pass)), alpha=0.8) + 
  labs(title="phasing pass vs width", 
       subtitle="City Mileage Grouped by Number of cylinders",
       caption="Source: mpg",
       x="Inv width",
       fill="phase pass")

g <- ggplot(tab, aes(log10(width)))
g + geom_density(aes(fill=factor(verdict)), alpha=0.8) + 
  labs(title="phasing pass vs width", 
       subtitle="City Mileage Grouped by Number of cylinders",
       caption="Source: mpg",
       x="Inv width",
       fill="phase pass")

g <- ggplot(tab, aes(verdict, log10(width)))
g + geom_boxplot()# + 
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5, 
               fill="red") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(title="Box plot + Dot plot", 
       subtitle="City Mileage vs Class: Each dot represents 1 row in source data",
       caption="Source: mpg",
       x="Class of Vehicle",
       y="City Mileage")


# variant allele frequency
cm = callmatrix
library(dplyr)


# vaf figure. ###########################################3
cm = callmatrix
#cm = callmatrixpass
cm$nhom = rowSums2(cm=='1|1')
cm$nhet  =rowSums2(cm=='1|0') + rowSums2(cm=='0|1')
cm$nref = rowSums2(cm=='0|0')
cm$ninv = (2*rowSums2(cm=='1|1')) + rowSums2(cm=='1|0') + rowSums2(cm=='0|1') + 
  (2*rowSums2(cm=='1|1_lowconf')) + rowSums2(cm=='1|0_lowconf') + rowSums2(cm=='0|1_lowconf')
cm$n = 2*((rowSums2(cm=='1|1')) + rowSums2(cm=='1|0') + rowSums2(cm=='0|1') + 
  (rowSums2(cm=='1|1_lowconf')) + rowSums2(cm=='1|0_lowconf') + rowSums2(cm=='0|1_lowconf') +
  (rowSums2(cm=='0|0_lowconf')) + (rowSums2(cm=='0|0')))

cm$ninv = cm$ninv + 


cm$af = cm$ninv / cm$n

cm2 = cm#[cm$n > 50,]
ggplot(data=cm2) + geom_histogram(aes(x=ninv))   # + scale_x_log10()

as.data.frame.matrix(table(cm$ninv))
aaa = as.data.frame(as.matrix(table(cm$ninv)))
aaa$n = as.numeric(row.names(aaa))




ggplot(data=aaa[aaa$n > 0,], aes(x=n, y=V1)) + geom_point() + geom_line()

ggplot(data=aaa[aaa$n > 0,], aes(x=n, y=V1)) + geom_point() + geom_line() + xlim(1,100) + scale_x_log10() +
  labs(title('Variant allele frequency'), x='Variant allele count', y='Inversion sites') + scale_y_log10()

ggplot(data=cm) + geom_histogram(aes(x=af)) 


######################################################3
#figure1
tab3 = tab#[tab$verdict == 'pass',]
companion_unique_samples = c('HG00514','HG00733','NA19240','HG02018','HG01573','GM19036')
# filter, filter
tab3a = left_join(tab3, cm[,c('chrom','start','end','simpleverdict')])
tab3b = tab3a[tab3a$simpleverdict %in% c('PASS','lowmap', 'AlwaysComplex'),]
tab3c = tab3b[!(tab3b$sample %in% companion_unique_samples),]
tab3 = tab3c

#tab3$GT = tab3$GTs
tab3$simpler = 'ccomplex'
tab3[tab3$GT %in% c('1|1', '1|1_lowconf'),]$simpler = 'aHOM'
tab3[tab3$GT %in% c('1|0', '1|0_lowconf','0|1', '0|1_lowconf'),]$simpler = 'bHET'
tab3[tab3$GT %in% c('0|0', '0|0_lowconf'),]$simpler = 'REF'
tab3[tab3$GT %in% c('noreads'),]$simpler = 'zeroreads'
tab3[tab3$ID %in% cm_detail_hom$ID,]$simpler = 'allMISO'
#tab3[(tab3$ID %in% cm_detail_hom$ID) & (tab3$sample == 'GM20509'),]$simpler = 'aaMISO'

#tab[endsWith(tab$GT, '_lowconf'),]$simpler = 'complex_lowconf'
#tab[tab$GT %in% c('0|0_lowconf', '1|0_lowconf', '0|1_lowconf', '1|1_lowconf'),]$simpler = 'simple_lowconf'
tab3 = tab3[tab3$simpler %in% c('aHOM','bHET','ccomplex', 'aaMISO','zeroreads'),]
#tab3[tab3$simpler %in%  c('aHOM','bHET','ccomplex'),]$simpler = 'Hom/Het/Complex'
data = as.data.frame.matrix(table(tab3$sample, tab3$simpler))
dat = melt(as.matrix(data))
library(ggbeeswarm)
g <- ggplot(data=dat, aes(x=X2, y=value)) +  geom_boxplot(aes(color=X2)) +
  # Plot first the fat dots
  geom_beeswarm(data=dat,
                size=2,
                cex=1.5) + 
  ylim(c(1,200)) + labs(x='Event class', y='Events per genome') +
  theme(axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold",
                                   size=14))
  
g

# inverted sequence per genome.
t4 = tab3 %>% group_by(sample,simpler) %>% mutate(nbases = sum(len)) %>% filter(row_number() == 1) 

t5 = t4[,c('sample','simpler','nbases')]
t6 = t5[t5$simpler %in% c('aHOM','bHET'),]
t6[t6$simpler == 'aHOM',]$nbases = 2 * t6[t6$simpler == 'aHOM',]$nbases
t7 = t6 %>% group_by(sample) %>% mutate(value = sum(nbases)) %>% filter(row_number() == 1) 
t7$x = 'sample'
t7 = as.data.frame(t7)
data = as.data.frame.matrix(table(t7$x, t7$value))
dat = melt(as.matrix(data))
t7$value = t7$value/1000 
g <- ggplot(data=t7, aes(x=x, y=value)) + 
  geom_boxplot(aes(color=x)) +
  geom_point(aes(x=x,y=value)) + 
  ylim(c(1,40)) +
  labs(x='', y='Inverted bases per diploid genome [Mb]') +
  theme(axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold",
                                   size=14))

g


#  accidental HW
library(HardyWeinberg)
CX = cm[,c('nref','nhet','nhom')]
HWTernaryPlot(CX,100,region=1,hwcurve=TRUE,vbounds=FALSE,vertex.cex=2)
  
dat$anc = 'UNK'
dat[dat$X1=='GM12329',]$anc = 'EUR'
  dat[dat$X1=='GM18534',]$anc = 'CHS'
  dat[dat$X1=='GM18939',]$anc = 'CHS'
  dat[dat$X1=='GM19036',]$anc = 'AFR'
  dat[dat$X1=='GM19650',]$anc = 'AMR'
  dat[dat$X1=='GM19983',]$anc = 'AFR'
  dat[dat$X1=='GM20509',]$anc = 'EUR'
  dat[dat$X1=='GM20847',]$anc = 'SAS'
  dat[dat$X1=='HG00096',]$anc = 'EUR'
  dat[dat$X1=='HG00171',]$anc = 'EUR'
  dat[dat$X1=='HG00512',]$anc = 'EAS'
  dat[dat$X1=='HG00513',]$anc = 'EAS'
  dat[dat$X1=='HG00514',]$anc = 'EAS'
  dat[dat$X1=='HG00731',]$anc = 'AMR'
  dat[dat$X1=='HG00732',]$anc = 'AMR'
  dat[dat$X1=='HG00733',]$anc = 'AMR'
  dat[dat$X1=='HG00864',]$anc = 'EAS'
  dat[dat$X1=='HG01114',]$anc = 'AMR'
  dat[dat$X1=='HG01505',]$anc = 'EUR'
  dat[dat$X1=='HG01573',]$anc = 'AMR'
  dat[dat$X1=='HG01596',]$anc = 'EAS'
  dat[dat$X1=='HG02011',]$anc = 'AFR'
  dat[dat$X1=='HG02018',]$anc = 'EAS'
  dat[dat$X1=='HG02492',]$anc = 'SAS'
  dat[dat$X1=='HG02587',]$anc = 'AFR'
  dat[dat$X1=='HG02818',]$anc = 'AFR'
  dat[dat$X1=='HG03009',]$anc = 'SAS'
  dat[dat$X1=='HG03065',]$anc = 'AFR'
  dat[dat$X1=='HG03371',]$anc = 'AFR'
  dat[dat$X1=='HG03683',]$anc = 'SAS'
  dat[dat$X1=='HG03732',]$anc = 'SAS'
  dat[dat$X1=='NA19238',]$anc = 'AFR'
  dat[dat$X1=='NA19239',]$anc = 'AFR'
  dat[dat$X1=='NA19240',]$anc = 'AFR'
                                 
#donut


donut = as.data.frame(as.matrix(table(tab3$simpler)))
donut$Prediction = row.names(donut)
donut = donut[c('REF','HET','HOM','complex','zeroreads'),]
donut$Prediction = c('aREF','bHET','cHOM','dcomplex','zeroreads')
row.names(donut) =  c('aREF','bHET','cHOM','dcomplex','zeroreads')
ypos1 = head(c(0,cumsum(rev(donut$V1))),-1)
ypos2 = cumsum(rev(donut$V1))
ypos = (ypos2 + ypos1)/2
  ggplot(donut, aes(x = 2, y = V1, fill = Prediction)) +
    geom_bar(stat = "identity", color = "white") +
    coord_polar(theta = "y", start = 0)+
    geom_text(aes(y = ypos, label = rev(V1)), color = "white")+
    theme_void()+
    xlim(0.5, 2.5)
  
# length plot
  
tab1 = tab[tab$sample=='HG00733',]

bins = 40

g = ggplot(data=tab1, aes(len)) + geom_histogram(aes(y=..density..), bins=bins, color = "black", fill = "grey") + scale_x_log10() + theme_minimal() + geom_density(aes(x=len), size=1)
g
g = ggplot(data=tab1, aes(len)) + geom_histogram(bins=bins, color = "black", fill = "grey") + xlim(c(0,1000)) + scale_x_log10() + geom_density(aes(x=len), size=1) 
g

tab4 = tab3
tab4[tab4$simpler == 'REF',]$simpler = 'eREF'
tab4[tab4$simpler == 'HET',]$simpler = 'dHET'
tab4[tab4$simpler == 'HOM',]$simpler = 'cHOM'
tab4[tab4$simpler == 'complex',]$simpler = 'bcomplex'
tab4[tab4$simpler == 'zeroreads',]$simpler = 'anoreads'

g = ggplot(data=tab4[tab4$sample=='HG00733',], aes(len)) + geom_histogram(bins=bins, aes(fill=simpler)) + scale_x_log10()  + theme_minimal() #+ geom_density(aes(x=len), size=1) 
g


ybreaks = seq(0,20,5) 
## On primary axis
g + scale_y_continuous("Counts", breaks = round(ybreaks / (bw * n_obs),3), labels = ybreaks)

## Or on secondary axis
g + scale_y_continuous("Density", sec.axis = sec_axis(
  trans = ~ . * bw * n_obs, name = "Counts", breaks = ybreaks))



ggplot(data=tab1) + geom_histogram(aes(x=len), bins=50) + scale_x_log10()
ggplot(data=tab1) + geom_histogram(aes(x=len), bins=60) + scale_x_log10()
