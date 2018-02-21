library(data.table)
library(assertthat)
library(pheatmap)
source("utils/sv_classifier_probs.R")
source("utils/sv_classifier_counts.R")



# counts = fread(paste("zcat","counts/simulation23-100000/100000_fixed.txt.gz"))
# info   = fread("counts/simulation23-100000/100000_fixed.info")
# strand = fread("strand_states/simulation23-100000/final.txt")
# segs   = fread("segmentation2/simulation23-100000/100000_fixed.few.txt") 

# counts = fread(paste("zcat","counts/simulation7-50000/50000_fixed.txt.gz"))
# info   = fread("counts/simulation7-50000/50000_fixed.info")
# strand = fread("strand_states/simulation7-50000/final.txt")
# segs   = fread("segmentation2/simulation7-50000/50000_fixed.few.txt")

counts = fread(paste("zcat",snakemake@input[["counts"]]))
info   = fread(snakemake@input[["info"]])
strand = fread(snakemake@input[["states"]])
segs   = fread(snakemake@input[["bp"]])






################################################################################
# Check input data 
#
# counts
assert_that("chrom" %in% colnames(counts),
            "start" %in% colnames(counts),
            "end"   %in% colnames(counts),
            "sample"%in% colnames(counts),
            "cell"  %in% colnames(counts))
counts <- counts[order(sample,cell,chrom,start,end),]
setkey(counts,sample,cell)

# info
assert_that("sample"%in% colnames(info),
            "cell"  %in% colnames(info),
            "nb_p"  %in% colnames(info),
            "nb_r"  %in% colnames(info),
            "nb_a"  %in% colnames(info),
            "pass1" %in% colnames(info))
info <- info[order(sample,cell),]
setkey(info,sample,cell)

# strand
assert_that("sample"%in% colnames(strand),
            "cell"  %in% colnames(strand),
            "chrom" %in% colnames(strand),
            "start" %in% colnames(strand),
            "end"   %in% colnames(strand),
            "class" %in% colnames(strand))
strand <- strand[order(sample,cell,chrom,start,end),]

# segs
assert_that("chrom" %in% colnames(segs),
            "bps"   %in% colnames(segs))
segs <- segs[order(chrom, bps),]

message("[SV classifier] Problem size: ", nrow(info), " cells x ", nrow(segs), " segments.")
################################################################################






# Kick out non-PASS cells     # To test, use something like  info[seq(1,nrow(info),3)]$pass1 = 0
if (nrow(info[pass1 != 1])> 0) message("[SV classifier] Kicking out ", nrow(info[pass1 != 1]), " low quality cells. ", nrow(info[pass1 == 1]), " remain.")
info <- info[pass1 == 1,]
counts <- counts[ paste(sample,cell) %in% info[,paste(sample,cell)] ]
assert_that(all(unique(counts[,.(sample,cell)]) == unique(info[,.(sample,cell)])))


# Get mean and median from count data
info <- merge(info, counts[, .(mean = mean((w+c)[class != "None"]), median = median((w+c)[class != "None"])), by = .(sample, cell)], by = c("sample","cell"))



# Prepare function to assess the strand_state of a cell/chrom/region
bins <- unique(counts[, .(chrom, start, end)])
get_strand_state <- function(sample_, cell_, chrom_, from_, to_) {
    x = strand[sample == sample_ & cell == cell_ & chrom == chrom_]
    assert_that(nrow(x)>0)
    if (nrow(x) == 1) return (x$class)
    min_pos = bins[chrom == chrom_]$start[from_]
    max_pos = bins[chrom == chrom_]$end[to_]
    x = x[start <= min_pos & end >= max_pos]
    if (nrow(x) == 1) return (x$class)
    return ("sce")
}


# Expand segments into {chrom, [from, to]}
segs[, from := shift(bps,fill = 0) + 1, by = chrom]
segs[, `:=`(to = bps, bps = NULL, k = NULL)]




# Do all the work:

message("[SV classifier] Preparing large table [segments + cells  x features]")
probs <- segs[,cbind(.SD,info[,.(sample,cell,nb_p,nb_r,medbin,mean)]), by = .(chrom,from)]

message("[SV classifier] Annotating strand-state (slow)")
probs[, state := get_strand_state(sample, cell, chrom, from, to), by = .(sample, cell, chrom, from, to)]

message("[SV classifier] Annotating expected coverage")
probs[, expected := (to - from +1)*mean, by = .(sample, cell, chrom, from, to)]

message("[SV classifier] Annotating observed W/C counts")
probs <- add_seg_counts(probs, counts)
probs[, scalar := 1]

message("[SV classifier] Annotating NB probabilities")
probs <- add_NB_probs(probs)

message("[SV classifier] Post-processing NB probabilities")
probs[state == "sce", `:=`(p_ref = 0, p_homInv = 0, p_hetInv = 0, p_hetDel = 0, p_homDel = 0)]
probs[,likelyhoodratio := p_ref - pmax(p_hetInv, p_hetDel, p_hetInv)]
probs[,obs_exp := (W+C)/expected]




# Model loci across all cells. 
# Each model allows only one type of SV in the locus
mod = probs[, data.table(model = c("ref","hetDel","homDel","hetInv","homInv","hetDup"),
                         loglik = c(sum(p_ref),
                                    sum(pmax(p_ref, p_hetDel)),
                                    sum(pmax(p_ref, p_homDel)),
                                    sum(pmax(p_ref, p_hetInv)),
                                    sum(pmax(p_ref, p_homInv)),
                                    sum(pmax(p_ref, p_hetDup))   ),
                         num    = c(sum(p_ref == pmax(p_ref, p_hetDel, p_homInv, p_hetInv, p_homDel, p_hetDup)),
                                    sum(p_hetDel == pmax(p_ref, p_hetDel)),
                                    sum(p_homDel == pmax(p_ref, p_homDel)),
                                    sum(p_hetInv == pmax(p_ref, p_hetInv)),
                                    sum(p_homInv == pmax(p_ref, p_homInv)),
                                    sum(p_hetDup == pmax(p_ref, p_hetDup))  )
                         )[order(loglik, decreasing = T)],
      by = .(chrom, from, to)]
MIN_CELLS = 2
mod = mod[, .SD[num>=MIN_CELLS][1,], by = .(chrom, from, to)]




# Apply the best model to the prob by overwriting probabilities
probs = merge(probs, mod[, .(chrom, from, to, model)], by = c("chrom","from","to"))
newprobs = probs

newprobs[model == "ref", p_ref := 0]
newprobs[model == "hetDup" & p_hetDup > p_ref, p_hetDup := 0]
newprobs[model == "hetDup" & p_hetDup <= p_ref, p_ref := 0]
newprobs[model == "hetDel" & p_hetDel > p_ref, p_hetDel := 0]
newprobs[model == "hetDel" & p_hetDel <= p_ref, p_ref := 0]
newprobs[model == "hetInv" & p_hetInv > p_ref, p_hetInv := 0]
newprobs[model == "hetInv" & p_hetInv <= p_ref, p_ref := 0]
newprobs[model == "homInv" & p_homInv > p_ref + 0.01, p_homInv := 0]
newprobs[model == "homInv" & p_homInv <= p_ref + 0.01, p_ref := 0]
newprobs[model == "homDel" & p_homDel > p_ref, p_homDel := 0]
newprobs[model == "homDel" & p_homDel <= p_ref, p_ref := 0]





# Before output, Rename columns
out = newprobs[, .(chrom,
                start = bins[from]$start, 
                end = bins[to]$end, 
                sample, 
                cell, 
                type = state, 
                w = W, 
                c = C, 
                p_ref     = exp(p_ref), 
                p_del_hom = exp(p_homDel), 
                p_del_h1  = exp(p_hetDel),
                p_del_h2  = exp(p_hetDel),
                p_inv_hom = exp(p_homInv),
                p_inv_h1  = exp(p_hetInv),
                p_inv_h2  = exp(p_hetInv),
                p_dup_h1  = exp(p_hetDup),
                p_dup_h2  = exp(p_hetDup)  )]
write.table(out, file = snakemake@output[[1]], row.names = F, col.names = T, quote=F, sep = "\t")

# heat map
#  X = dcast(probs[grepl('^chr3', chrom)], chrom + from + to ~ sample + cell, value.var = "likelyhoodratio")
#  Y = as.matrix(X[, 4:ncol(X), with = F])
#  assert_that(all(!is.infinite(Y)))
#  rownames(Y) = paste(X$chrom, X$from, X$to, sep = "_")
#  Y[Y > 100] = 100
#  Y[Y < -100] = -100
#  pheatmap(Y, breaks = c(-101,-5,0,5,101), color = c("blue","turquoise","orange","red"))

