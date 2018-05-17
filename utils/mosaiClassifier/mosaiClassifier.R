library(dplyr)
library(data.table)
library(assertthat)
source("utils/mosaiClassifier/getStrandStates.R")
source("utils/mosaiClassifier/getCountsPerSegment.R")
source("utils/mosaiClassifier/generateHaploStates.R")
source("utils/mosaiClassifier/getDispParAndSegType.R")
source("utils/mosaiClassifier/haploAndGenoName.R")


dir <- "/home/maryam/research/hackathons/troubleshooting/simulatedData/simulation5-100000/"
binRCfile <- paste0(dir, "100000_fixed.txt.gz")
BRfile <- paste0(dir, "100000_fixed.few.txt")
infoFile <- paste0(dir, "100000_fixed.info")
stateFile <- paste0(dir, "final.txt")
counts <- fread(paste("zcat", binRCfile))
info <- fread(infoFile)
strand <- fread(stateFile)
segs <- fread(BRfile)

mosaiClassifierPrepare <- function(counts, info, strand, segs) {

  ##############################################################################
  # Check input data
  #
  assert_that(is.data.table(counts),
              "chrom" %in% colnames(counts),
              "start" %in% colnames(counts),
              "end"   %in% colnames(counts),
              "sample"%in% colnames(counts),
              "cell"  %in% colnames(counts)) %>% invisible
  setkey(counts, sample, cell, chrom, start, end)

  assert_that(is.data.table(info),
              "sample"%in% colnames(info),
              "cell"  %in% colnames(info),
              "nb_p"  %in% colnames(info),
              "nb_r"  %in% colnames(info),
              "nb_a"  %in% colnames(info),
              "pass1" %in% colnames(info)) %>% invisible
  setkey(info,sample,cell)

  assert_that(is.data.table(strand),
              "sample"%in% colnames(strand),
              "cell"  %in% colnames(strand),
              "chrom" %in% colnames(strand),
              "start" %in% colnames(strand),
              "end"   %in% colnames(strand),
              "class" %in% colnames(strand)) %>% invisible
  setkey(strand, sample, cell, chrom, start, end)

  # segs
  assert_that(is.data.table(segs),
              "chrom" %in% colnames(segs),
              "bps"   %in% colnames(segs)) %>% invisible
  setkey(segs, chrom, bps)


  # Kick out non-PASS cells
  if (nrow(info[pass1 != 1])> 0) {
    message("[SV classifier] Kicking out ",
            nrow(info[pass1 != 1]),
            " low quality cells. ",
            nrow(info[pass1 == 1]),
            " remain.")
    info <- info[pass1 == 1,]
    counts <- counts[ paste(sample,cell) %in% info[,paste(sample,cell)] ]
  }
  # Check that set of cells in counts is the same as in info
  assert_that(all(unique(counts[,.(sample,cell)]) == unique(info[,.(sample,cell)]))) %>% invisible



  ##############################################################################
  # Estimation mean read count per bin for dispersion parameters
  #
  message("[MosaiClassifier] Problem size: ",
          nrow(info),
          " cells x ",
          nrow(segs),
          " segments.")

  # Get trimmed mean of counts per bin to estimate the r parameter
  # When calculating this, ignore the "None" bins
  info <- suppressWarnings(
          merge(info,
                counts[, .(mean = mean((w+c)[class != "None"], trim = 0.05)),
                       by = .(sample, cell)],
                by = c("sample","cell")) )


  ##############################################################################
  # Expand table to (all cells) x (all segments)
  #
  # add a "from" column which contians the "to" breakpoint from the prev. segment each
  segs[, from := shift(bps,fill = 0) + 1, by = chrom]

  # rename the "bps" column to "to"
  segs[, `:=`(to = bps, bps = NULL, k = NULL)]

  # Add coordinates
  addPositions(segs, counts)

  # Expand the table to (cells) x (segments) and annotate with NB params
  # --> take each row in "segs" and cbind it to a whole "info" table
  probs <- segs[,
                cbind(.SD, info[,.(sample, cell, nb_p, mean)]),
                by = .(chrom,from)]


  ##############################################################################
  # Annotate each segment and cell with the strand state
  #
  message("[MosaiClassifier] Annotating strand-state")
  probs = addStrandStates(probs, strand)

  ##############################################################################
  # Annotate the observed and expected counts in each segment / cell
  #
  message("[MosaiClassifier] Annotating expected coverage")
  probs[, expected := (to - from +1)*mean, by = .(sample, cell, chrom, from, to)]

  message("[MosaiClassifier] Annotating observed W/C counts")
  probs <- addCountsPerSegment(probs, counts)
  probs[, scalar := 1]

  return(probs)
}



mosaiClassifierCalcProbs <- function(probs, maximumCN=4, haplotypeMode=F, alpha=0.05) {

  assert_that(is.data.table(probs))
  # check the colnames
  # defining the vector of all possible haplotypes
  hapStatus <- NULL
  for (j in 0:maximumCN)
  {
    hapStatus <- c(hapStatus, allStatus(3,j))
  }
  for (j in 1:length(hapStatus))
  {
    hapStatus[j] <- paste(decodeStatus(hapStatus[j]), collapse = '')
  }
  
  # creating a datatable containing all possible combinations of strand states and haplotypes,
  # and setting their segTypes
  hapStrandStates <- data.table()
  for (st in c("CC","WW","WC","CW"))
  {
    hapStrandStates <- rbindlist(list(hapStrandStates, 
                                      data.table(class=st, haplotype=hapStatus, 
                                                 segtype=t(sapply(hapStatus, function(x) getSegType(st, x))))))
  }
  # naming third and forth columns
  colnames(hapStrandStates)[3:4]=c("Wcn", "Ccn")
  # adding haplotype and genotype name columns
  hapStrandStates[,haplo_name:=.(sapply(haplotype, get_hap_name))]
  hapStrandStates[,geno_name:=.(sapply(haplo_name, haplo_to_geno_name))]
  # sort based on state
  setkey(hapStrandStates, class)
  
  ##### mering probs and haplotype strand states
  # kick out the segs with sces
  probs <- probs[class!="?"]
  
  probs <- merge(probs, 
                hapStrandStates,
                by = "class",
                allow.cartesian = T)
  ###########
  # compute dispersion parameter (nb_r) column
  probs[, nb_r:=.(nb_p*mean/(1-nb_p))]
  # reshuffling the columns
  probs <- probs[,.(sample, cell, chrom, start, end, from, to, nb_p, nb_r, mean, class, expected,
           W, C, scalar, haplotype, Wcn, Ccn, haplo_name, geno_name)]
  
  # computing dispersion parameters seperately for each segment and W and C counts
  probs <- add_dispPar(probs, alpha)
  
  # compute NB haplotype likelihoods
  probs[, nb_hap_ll := dnbinom(W, size = disp_w, prob = nb_p)
        *dnbinom(C, size = disp_c, prob = nb_p)]
  
  # computing sister haplotype (haplotype with the same genotype) for each haplotype
  sister.haps <- sapply(hapStatus, sisterHaplotype)
  sister.hap.pos <- match(sister.haps, hapStatus)
  # compute the set of symmetric haplotypes (haplotypes that are equal to their sister haplotype)
  symmetric.haps <- hapStatus[which(sister.haps==hapStatus)]
  
  # averaging the nb probs of sister haplotypes, when haplotype specific strand states are not known
  # adding genotype likelihoods, if haplotype mode is false
  if (!haplotypeMode)
  {
    probs[,nb_hap_ll:=.((nb_hap_ll+nb_hap_ll[sister.hap.pos])/2), by=.(sample, cell, chrom, from, to)]
    
    # computing genotype likelihoods
    probs[,nb_gt_ll:=.(nb_hap_ll+nb_hap_ll[sister.hap.pos]), by=.(sample, cell, chrom, from, to)]
    # deviding the gt likelihoods of symmetric haplotypes by 2
    probs[haplotype %in% symmetric.haps, nb_gt_ll:=.(nb_gt_ll/2)]
  }
  
  # TODO export the prob table to some output file
  
  return(probs)
}

mosaiClassifierPostProcessing <- function(probs, haplotypeMode=F, regularizationFactor=1e-10)
{
  assert_that(is.data.table(probs))
  # check the colnames
  
  # testing if there are some segments with zero probability for all haplotypes
  segs_max_hap_nb_probs <- probs[,
                                 .(sample, chrom, cell, from, to, max_nb_hap_ll=rep(max(nb_hap_ll), .N)), 
                                 by=.(sample, chrom, cell, from, to)]
  message(paste("the number of segments with 0 prob for all haplotypes = ", 
                segs_max_hap_nb_probs[max_nb_hap_ll==0, .N]))
  
  # add prior probs to the table
  probs[,prior:=100L]
  probs[haplo_name=="ref_hom",prior:=200L]
  probs[haplo_name=="complex",prior:=1L]
  
  # compute the posteriori probs (add new columns)
  probs[,nb_hap_pp:=.(nb_hap_ll*prior)][,nb_gt_pp:=.(nb_gt_ll*prior)]
  
  # set a uniform prob on sce segs and the segs_max_hap_nb_probs=0
  probs[segs_max_hap_nb_probs$max_nb_hap_ll==0,nb_hap_pp:=1L]
  probs[class=="?", nb_hap_pp:=1L]
  
  # normalizing nb_hap_pp and nb_gt_pp to 1 per sample, cell, and segment
  probs[, nb_hap_pp := nb_hap_pp/sum(nb_hap_pp), by=.(sample, cell, chrom, from, to)]
  probs[, nb_gt_pp := nb_gt_pp/sum(nb_gt_pp), by=.(sample, cell, chrom, from, to)]
  
  # regularizing nb_hap_ll to set the min possible likelihood to a constant small number
  probs[,nb_hap_pp:=.((regularizationFactor/length(hapStatus))+nb_hap_pp*(1-regularizationFactor))]
  probs[,nb_gt_pp:=.((regularizationFactor/length(hapStatus))+nb_gt_pp*(1-regularizationFactor))]
  
  # converting to simple haplotype prob table
  simp.probs <- probs[haplo_name!="complex"]
  
  # normalizing hap and gt probs in simp.probs
  probs[, nb_hap_pp := nb_hap_pp/sum(nb_hap_pp), by=.(sample, cell, chrom, from, to)]
  probs[, nb_gt_pp := nb_hap_pp/sum(nb_gt_pp), by=.(sample, cell, chrom, from, to)]
  
  # dcasting: converting the table from long to wide format based on the haplotype names
  
  # calling SVs (It should be included in Sascha's code)
}
