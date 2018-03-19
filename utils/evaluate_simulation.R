library(data.table)
library(assertthat)
library(ggplot2)
library(cowplot)
library(scales)

MIN_LLR = 0.5   # min. log likelihood ratio of SV prob. over REF prob.
MIN_RCO = 0.7   # min. reciprocal overlap of true and detected SV

#sv_input = "sv_probabilities/simulation7-100000/100000_fixed.few/probabilities.txt"
#real_input = "simulation/variants/genome7-100000.txt"
sv_input = snakemake@input[["prob"]]
real_input = snakemake@input[["simul"]]


svs = fread(sv_input)
real = fread(real_input)


### call SVs from SV probabilities

# rename p_ref column to be distinguished from SV classes
svs[, `:=`(ref_prob = p_ref, p_ref = NULL) ]

# melt to one entry per SV class
svs = melt.data.table(svs, measure.vars = colnames(svs)[grepl('^p_', colnames(svs))], variable.name = "SV_class", value.name = "SV_prob", variable.factor = FALSE)

# Rename called SV classes
unique(svs$SV_class)
svs <- svs[SV_class %in% c("p_inv_hom", "p_inv_h1", "p_inv_h2", "p_del_hom", "p_del_h1",  "p_del_h2", "p_dup_hom", "p_dup_h1", "p_dup_h2", "p_idup_h1", "p_idup_h2") ]
svs[SV_class == "p_inv_hom", SV_class := "hom_inv"]
svs[SV_class == "p_inv_h1", SV_class := "het_inv"]
svs[SV_class == "p_inv_h2", SV_class := "het_inv"]
svs[SV_class == "p_del_hom", SV_class := "hom_del"]
svs[SV_class == "p_del_h1", SV_class := "het_del"]
svs[SV_class == "p_del_h2", SV_class := "het_del"]
svs[SV_class == "p_dup_hom", SV_class := "hom_dup"]
svs[SV_class == "p_dup_h1", SV_class := "het_dup"]
svs[SV_class == "p_dup_h2", SV_class := "het_dup"]
svs[SV_class == "p_idup_h1", SV_class := "inv_dup"]
svs[SV_class == "p_idup_h2", SV_class := "inv_dup"]


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

# Merge all combinations of loci by chromosoems ...
combined = merge(LOCI[, .(chrom, start, end)],
                 called_loci,
                 by = "chrom",
                 all =T,
                 suffixes = c("",".call"),
                 allow.cartesian = T)

# ... to then filter those that have a reciprocal overlap of at least MIN_RCO
combined = combined[!is.na(start) & end.call > start & start.call < end]
combined[, recovl := (pmin(end, end.call) - pmax(start, start.call)) / (pmax(end, end.call) - pmin(start, start.call))]

# annotate simulated loci with matched real loci, if they exist.
LOCI = merge(LOCI,
             combined[recovl >= MIN_RCO],
             by = c("chrom","start","end"), all.x = T)
LOCI[, id := 1:.N]
assert_that(all(LOCI[,.(chrom, start, end)] == unique(real[, .(chrom, start, end)])))


### Determine sensitivity (recall)
# annotate "real" and "svs" with locus ids for easier overlap
# Note: svs contains also non-SV loci. Always filter by SV_llr >= MIN_LLR
#
real = merge(real, LOCI[, .(chrom, start, end, id)], by =c("chrom","start","end"), all = T)
svs  = merge(svs, LOCI[, .(chrom, start = start.call, end = end.call, id)], by =c("chrom","start","end"), all.x = T)




### EVALUATION PART

# recall    = svs correctly detected / svs simulated
#
recall = merge(real[, .(id, sample, cell, simulated = SV_type)], 
               svs[!is.na(id) & SV_llr >= MIN_LLR, .(id, sample, cell, called = SV_class)], 
               by = c("id","sample","cell"), all.x = T)
recall = recall[, .(recall.correct = sum(simulated[!is.na(called)] == called[!is.na(called)]), 
                    recall.any = sum(!is.na(called)), 
                    recall.total = .N), 
                by = .(id)]
LOCI = merge(LOCI, recall, by = "id")



# Precision
# precision = svs correctly detected / all detected svs
# First get all calls that do not match simulated SV
spurious_calls = svs[SV_llr >= MIN_LLR,
                     .N,
                     by = .(chrom, start, end, id, SV_class)] # all loci, all classes
# There are also some calls that overlap SVs (by at least 0.1 %) without fully matching them
twighlight_loci = combined[recovl >= 0.001,
                           .(chrom, start = start.call, end = end.call, twighlight = T)]
# add these overlapping loci to the list of spurious calls
spurious_calls = merge(spurious_calls,
                       twighlight_loci,
                       by = c("chrom","start","end"),
                       all.x = T)
message("Summary of spurious loci either matching or at least overlapping simulated SVs")
table(factor(is.na(spurious_calls$id), c(T,F), c("No SV match", "SV match")),
      factor(is.na(spurious_calls$twighlight), c(T,F), c("no SV overlap", "SV overlap")))
assert_that(nrow(spurious_calls[!is.na(id) & is.na(twighlight),]) == 0)
# Only keep the actual spurious calls, i.e. the ones that do not match a simulated SV
# (because they have no ID) OR the ones that are in the twighlight zone (!i.sna(twighlight))
spurious_calls = spurious_calls[is.na(twighlight) | is.na(id),]
spurious_calls[, category := ifelse(is.na(twighlight), "really off", "some overlap with true SVs")]
spurious_calls[, twighlight := NULL]
spurious_calls[, size := factor((end - start) < 1e6, c(T,F), c("< 1Mb", ">= 1Mb"))]





### PLOTS ###



# Plot: number of breakpoints that were detected
format_Mb   <- function(x) {paste(round(x/1e6,1), "Mb")}
plt_obj = LOCI[, .(size = end-start, VAF, SV_type, detected = !is.na(start.call))]
plt1 <- ggplot(plt_obj) + 
  geom_point(size = 3, aes(size, VAF, col = detected)) +
  xlab("SV size") +
  scale_color_discrete(name = paste("based on",MIN_RCO*100,"% rec. ovl")) +
  theme(legend.position = "bottom") +
  scale_x_log10(labels = format_Mb, breaks = c(1e4,3e4,1e5,3e5,1e6,3e6,1e7)) +
  annotation_logticks(side = "b") +
  ggtitle(paste0("Were the breakpoints detected? ", nrow(plt_obj[detected==TRUE]), "/", nrow(plt_obj), " based on ", MIN_RCO*100, "% rec. ovlerap"))



# Plot: Recall on detected loci
plt_obj = LOCI[!is.na(start.call), .(size = end - start, VAF, SV_type, n = recall.total, k = recall.correct)]
plt2 = ggplot(plt_obj) + 
  geom_point(aes(size, VAF, col = k/n), size = 3) +
  scale_color_gradientn(name = "Recall (correct SV type)",
                        values = c(0,0.1,0.9,1),
                        colours = c("black","dodgerblue","darkorange","red")) +
  facet_wrap(~SV_type, nrow = 3) +
  theme(legend.position = "bottom", legend.key.width = unit(2,"cm")) +
  xlab("SV size") +
  scale_x_log10(labels = format_Mb, breaks = c(1e4,1e5,1e6,1e7)) +
  annotation_logticks(side = "b") +
  ggtitle("Recall")



# 3: Plot Number of spurious calls per locus
plt3 = ggplot(spurious_calls) +
  geom_histogram(aes(N, ), binwidth = 1) +
  facet_wrap(~ category + size, ncol = 1) +
  ggtitle(paste0("Spurious SV calls(", nrow(unique(spurious_calls[,.(chrom,start,end)])), " loci, ", sum(spurious_calls$N), " calls)")) +
  xlab("Num. cells with SV prediction at locus")



# 4: SV classes of spurious calls
plt4 = ggplot(spurious_calls) +
  geom_bar(aes(SV_class)) +
  facet_wrap(~ category + size, ncol = 1) +
  ggtitle(paste0("Wrongly called SV classes"))



final = ggdraw() +
  draw_label(sv_input, x=0.5, y=.97, hjust = 0) +
  draw_plot(plt1,     0, 0.7, 0.5, 0.3) +
  draw_plot(plt2,     0,   0, 0.5, 0.7) +  # recall
  draw_plot(plt3,   0.5,   0, 0.25, 0.7) +
  draw_plot(plt4,  0.75,   0, 0.25, 0.7)

ggsave(snakemake@output[[1]], final, width = 29.7/1.3, height = 21/1.3)