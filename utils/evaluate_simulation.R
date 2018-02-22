library(data.table)
library(assertthat)
library(ggplot2)
library(cowplot)

MIN_LLR = 0.5   # min. log likelihood ratio of SV prob. over REF prob.
MIN_RCO = 0.8   # min. reciprocal overlap of true and detected SV

sv_input = "sv_probabilities/simulation7-100000/100000_fixed.few/probabilities.txt"
real_input = "simulation/variants/genome7-100000.txt"
sv_input = snakemake@input[["prob"]]
real_input = snakemake@input[["simul"]]




svs = fread(sv_input)
real = fread(real_input)


### call SVs from SV probabilities

# rename p_ref column to be distinguished from SV classes
svs[, `:=`(ref_prob = p_ref, p_ref = NULL) ]

# melt to one entry per SV class
svs = melt.data.table(svs, measure.vars = colnames(svs)[grepl('^p_', colnames(svs))], variable.name = "SV_class", value.name = "SV_prob", variable.factor = FALSE)

# cast again to one entry per locus and cell, now with 
setkey(svs, chrom, start, end, sample, cell)
svs = svs[, .(SV_class = SV_class[which.max(SV_prob)],
      SV_llr   = log(max(SV_prob)) - log(ref_prob[1])), 
  by = .(chrom, start, end, sample, cell)]



### Find reciprocal overlap of SV calls and true SVs

# list of simulated loci (only chrom, start, end)
LOCI = real[, .(VAF = .N), by = .(chrom, start, end, SV_type)]

# list of predicted loci (chrom, start, end) --> only use real SV calls here (MIN_LLR criterion) !
called_loci = unique(svs[SV_llr >= MIN_LLR, .(chrom, start, end)])

# Merge all combinations of loci ...
combined = merge(LOCI[, .(chrom, start, end)], called_loci, by = "chrom", all =T, suffixes = c("",".call"), allow.cartesian = T)

# ... to then filter those that have a reciprocal overlap of at least MIN_RCO
combined = combined[!is.na(start) & end.call > start & start.call < end & (pmin(end, end.call) - pmax(start, start.call)) / (pmax(end, end.call) - pmin(start, start.call)) >= MIN_RCO]

# annotate simulated loci with matched real loci, if they exist.
LOCI = merge(LOCI, combined, by = c("chrom","start","end"), all.x = T)
LOCI[, id := 1:.N]
assert_that(all(LOCI[,.(chrom, start, end)] == unique(real[, .(chrom, start, end)])))


# Rename called SV classes
unique(svs$SV_class)
svs[SV_class == "p_inv_hom", SV_class := "hom_inv"]
svs[SV_class == "p_inv_h1", SV_class := "het_inv"]
svs[SV_class == "p_inv_h2", SV_class := "het_inv"]
svs[SV_class == "p_del_hom", SV_class := "hom_del"]
svs[SV_class == "p_del_h1", SV_class := "het_del"]
svs[SV_class == "p_del_h2", SV_class := "het_del"]
svs[SV_class == "p_dup_hom", SV_class := "hom_dup"]
svs[SV_class == "p_dup_h1", SV_class := "het_dup"]
svs[SV_class == "p_dup_h2", SV_class := "het_dup"]



### Determine sensitivity (recall)
# annotate "real" and "svs" with locus ids for easier overlap
real = merge(real, LOCI[, .(chrom, start, end, id)], by =c("chrom","start","end"), all = T)
svs  = merge(svs, LOCI[, .(chrom, start = start.call, end = end.call, id)], by =c("chrom","start","end"), all.x = T)



# recall    = svs correctly detected / svs simulated
# precision = svs correctly detected / all detected svs 
#
recall = merge(real[, .(id, sample, cell, simulated = SV_type)], 
               svs[!is.na(id) & SV_llr >= MIN_LLR, .(id, sample, cell, called = SV_class)], 
               by = c("id","sample","cell"), all.x = T)
recall = recall[, .(recall.correct = sum(simulated[!is.na(called)] == called[!is.na(called)]), 
                    recall.any = sum(!is.na(called)), 
                    recall.total = .N), 
                by = .(id)]
LOCI = merge(LOCI, recall, by = "id")


prec = svs[SV_llr >= MIN_LLR, .N, by = .(chrom, start, end, id, SV_class)]




# Plot: number of breakpoints that were detected
plt_obj = LOCI[, .(size = end-start, VAF, SV_type, detected = !is.na(start.call))]
plt1 <- ggplot(plt_obj) + 
  geom_point(size = 2, aes(size, VAF, col = detected)) + 
  xlab("SV size") + 
  scale_color_discrete(name = paste("based on",MIN_RCO*100,"% rec. ovl")) + 
  theme(legend.position = "bottom") +
  ggtitle(paste0("Were the breakpoints detected? ", nrow(plt_obj[detected==TRUE]), "/", nrow(plt_obj)))


# Plot: Recall on detected loci
plt_obj = LOCI[!is.na(start.call), .(size = end - start, VAF, SV_type, n = recall.total, k = recall.correct)]
plt2 = ggplot(plt_obj) + 
  geom_point(aes(size, VAF, col = k/n), size = 1.5) +
  scale_color_gradientn(name = "Recall (correct SV type)",
                        values = c(0,0.1,0.9,1),
                        colours = c("black","dodgerblue","darkorange","red")) +
  facet_wrap(~SV_type, nrow = 1) +
  theme(legend.position = "bottom", legend.key.width = unit(2,"cm")) +
  xlab("SV size") +
  ggtitle("Recall")
plt2


plt_obj = prec[is.na(id), .N, by = .(chrom,start,end)]
plt_obj[, size := factor((end - start) < 1e6, c(T,F), c("< 1Mb", ">= 1Mb"))]
plt3 = ggplot(plt_obj) + 
  geom_histogram(aes(N), binwidth = 1) + 
  facet_wrap(~size) +
  ggtitle(paste0("SV calls outside simulated region (", sum(plt_obj$N), ")")) + 
  xlab("Num. calls at locus")
plt3

plt_obj = prec[is.na(id)]
plt4 = ggplot(plt_obj) + 
  geom_bar(aes(SV_class)) + 
  ggtitle(paste0("Wrongly called SV classes"))
plt4


final = ggdraw() + 
  draw_label(sv_input, x=0.5, y=.97) + 
  draw_plot(plt1,     0, 0.4, 0.33, 0.4) + 
  draw_plot(plt2,     0,   0,    1, 0.4) +
  draw_plot(plt3,  0.33, 0.4, 0.33, 0.4) +
  draw_plot(plt4,  0.66, 0.4, 0.33, 0.4)
final
ggsave(snakemake@output[[1]], final, width = 29.7/1.5, height = 21/1.5)
