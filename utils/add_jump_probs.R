library(GenomicRanges)
library(ggplot2)
#library(rgl) <- remove, because it caused problems for Docker/Singularity

dir <- "/local/home/maryam/research/hackathons/simulation-29May-2018/pipeline/"
binRCfile <- paste0(dir, "counts/simulation0-50000/50000_fixed.txt.gz")
BRfile <- paste0(dir, "segmentation2/simulation0-50000/50000_fixed.many.txt")
infoFile <- paste0(dir, "counts/simulation0-50000/50000_fixed.info")
stateFile <- paste0(dir, "strand_states/simulation0-50000/final.txt")
trueSVfile <- paste0(dir, "simulation/variants/genome0-50000.txt")
counts <- fread(paste("zcat", binRCfile))
info <- fread(infoFile)
strand <- fread(stateFile)
segs <- fread(BRfile)
trueSVs <- fread(trueSVfile)

probs <- mosaiClassifierPrepare(counts, info, strand, segs)
probs <- mosaiClassifierCalcProbs(probs)
probs <- mosaiClassifierPostProcessing(probs, haplotypeMode = T)
probs <- forceBiallelic(probs)
SV.calls <- makeSVCallSimple(probs, full.calls=T)

# adding genotype classes
probs[, geno_class:=haplo_code_to_geno_class(haplotype)]


getJumpProb <- function(probs, aggregateCells=T)
{
  # computing geno_class probs
  jump.probs <- probs[,.(class=class[1], nb_p=nb_p[1], expected=expected[1],
                         allele=allele[1], geno_class_pp = sum(nb_hap_pp)), 
                      by = .(sample, cell, chrom, start, end, geno_class)]
  
  # computing jump probs
  # sorting the probs table
  setkey(jump.probs, sample, chrom, start, end, cell, geno_class)
  # shift the probs column by size of #cells x #geno_classes(4)
  numCells <- length(unique(probs$cell))
  numGenoClasses <- 4
  
  # output jump probs
  jump.probs <- jump.probs[, next_seg_geno_class_pp:=
                             data.table::shift(geno_class_pp, n = numCells*numGenoClasses
                                               , fill=NA, type = "lead"),
                           by = chrom]
  # computing jump probs per segment per cell
  jump.probs <- jump.probs[,
                           .(class=class[1], nb_p=nb_p[1], expected=expected[1],
                             allele=allele[1], jump_pp=1-sum(geno_class_pp*next_seg_geno_class_pp)),
                           by = .(sample, chrom, start, end, cell)]
  
  # compute jump probs per segment
  jump.probs.segs <- jump.probs[, .(class=class[1], nb_p=nb_p[1], expected=expected[1],
                                    allele=allele[1], log_jump_pp=sum(log(jump_pp))),
                                by = .(sample, chrom, start, end)]
  
  # remove NA values
  jump.probs.segs <- jump.probs.segs[!is.na(log_jump_pp)]
  
  if (aggregateCells){
    return(jump.probs)
  } else {
    return(jump.probs.segs)
  }
}

# this function adds a column to the SV.calls table showing the number of similar 
# SV calls among the cells for each segment and its next segment
getSVconsistency <- function(SV.calls)
{
  setkey(SV.calls, sample, chrom, start, end, cell)
  numCells <- length(unique(probs$cell))
  
  # adding the column for SV calls of the next segments
  SV.calls[, next_seg_sv_call_haplotype:=
             data.table::shift(sv_call_haplotype, n = numCells, 
                             fill=NA, type = "lead"),
           by = chrom]
  
  SV.calls[!is.na(next_seg_sv_call_haplotype), 
           consistantCells:=length(which(sv_call_haplotype==next_seg_sv_call_haplotype)),
         by=.(sample, chrom, start, end)]
  
  SV.calls.consistency <- SV.calls[, .(cell=cell[1], class=class[1], 
                                     consistantCells=consistantCells[1]),
         by=.(sample, chrom, start, end)]
  
  return(SV.calls.consistency[!is.na(consistantCells)])
}

getAggProbDiff <- function(probs)
{
  setkey(probs, sample, chrom, start, end, cell)
  numCells <- length(unique(probs$cell))
  numHaps <- length(unique(probs$haplotype))
  
  # adding the column for agg hap prob of the next segments
  probs[, next_seg_agg_hap_pp:=
          data.table::shift(agg_hap_pp, n = numCells*numHaps, fill=NA, type = "lead"),
        by = chrom]
  
  probs <- probs[!is.na(next_seg_agg_hap_pp), 
                 .(agg_hap_log_pp_diff=sum(exp(agg_hap_pp+next_seg_agg_hap_pp))),
                 by=.(sample, chrom, start, end)]
  
  return(probs)
}

getTrueSV_BRs <- function(probs, trueSVs)
{
  # adding vaf column to trueSVs
  trueSVs[, vaf:=.N, by=.(chrom, start, end)]
  
  # labeling each segment based on whether the end breakpoint is true or false
  # creating a GRanges object for trueSVs and predicted SVs
  true.SVs.gr <- GRanges(unique(trueSVs[,.(chrom, start, end, vaf)]))
  
  segs.gr <- GRanges(unique(probs[,.(chrom, start, end)]))
  
  # finding overlaps between the true and the detected SVs
  ovp <- findOverlaps(segs.gr, true.SVs.gr)
  
  # create a data.table of the overlapping segments
  overlap <- data.table(chrom=as.character(seqnames((segs.gr[queryHits(ovp)]))),
                        pred.end=end(segs.gr[queryHits(ovp)]), 
                        true.start=start(true.SVs.gr[subjectHits(ovp)]), 
                        true.end=end(true.SVs.gr[subjectHits(ovp)]),
                        vaf=true.SVs.gr$vaf[subjectHits(ovp)])
  
  # label a predicted end as false if it doesn't match (distance less than 0.5*binSize) any start and end in the true SVs
  bin.size <- median(counts[,end-start])
  
  overlap[, trueBR:=min(abs(pred.end-true.start), abs(pred.end-true.end))<=bin.size/2
          , by =.(pred.end, true.start, true.end)]
  
  # remove repetitive predicted segments
  overlap <- overlap[,
                     {if (all(!trueBR)) {head(.SD,1)} else { head(.SD[trueBR],1) }},
                     by=.(chrom, pred.end)]
  return(overlap)
}

jump.probs.segs <- getJumpProb(probs)
trueBRs <- getTrueSV_BRs(probs, trueSVs)
SV.calls.consistency <- getSVconsistency(SV.calls)
aggP.consist <- getAggProbDiff(probs)

names(trueBRs)[2] <- "end"
# merge this trueBR column with the jump probs table
jump.probs.segs <- merge(trueBRs[,.(chrom, end, trueBR, vaf)], jump.probs.segs,
                         all.y = T, 
                         allow.cartesian=T,
                   by = c("chrom", "end"))[is.na(trueBR), trueBR:=F][is.na(vaf), vaf:=0]

SV.calls.consistency <- merge(trueBRs[,.(chrom, end, trueBR, vaf)], SV.calls.consistency,
                        all.y = T, 
                        allow.cartesian=T,
                        by = c("chrom", "end"))[is.na(trueBR), trueBR:=F][is.na(vaf), vaf:=0]

aggP.consist <- merge(trueBRs[,.(chrom, end, trueBR, vaf)], aggP.consist,
                              all.y = T, 
                              allow.cartesian=T,
                              by = c("chrom", "end"))[is.na(trueBR), trueBR:=F][is.na(vaf), vaf:=0]


# plotting
ggplot(data=jump.probs.segs, aes(x=log_jump_pp, y=end-start, col=trueBR)) + geom_point() + facet_grid(~trueBR)
ggplot(data=jump.probs.segs, aes(x=log_jump_pp, y=vaf, col=trueBR)) + geom_point() + facet_grid(~trueBR)

ggplot(data=SV.calls.consistency, aes(x=consistantCells, y=end-start, col=trueBR)) + geom_point()
ggplot(data=SV.calls.consistency, aes(x=consistantCells, y=vaf, col=trueBR)) + geom_point()

ggplot(data=aggP.consist, aes(x=agg_hap_log_pp_diff, y=end-start, col=trueBR)) + geom_point() + facet_grid(~trueBR)
ggplot(data=aggP.consist, aes(x=agg_hap_log_pp_diff, y=vaf, col=trueBR)) + geom_point()
