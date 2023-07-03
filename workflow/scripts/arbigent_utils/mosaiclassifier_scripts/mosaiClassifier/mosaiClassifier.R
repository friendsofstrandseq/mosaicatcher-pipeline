suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
source("workflow/scripts/arbigent_utils/mosaiclassifier_scripts/mosaiClassifier/getStrandStates.R")
source("workflow/scripts/arbigent_utils/mosaiclassifier_scripts/mosaiClassifier/getCountsPerSegment.R")
source("workflow/scripts/arbigent_utils/mosaiclassifier_scripts/mosaiClassifier/generateHaploStates.R")
source("workflow/scripts/arbigent_utils/mosaiclassifier_scripts/mosaiClassifier/getDispParAndSegType.R")
source("workflow/scripts/arbigent_utils/mosaiclassifier_scripts/mosaiClassifier/haploAndGenoName.R")

suppressMessages(library(GenomicRanges)) 
# This is to test whether the given strand states
# are in fact disjoint intervals !


mosaiClassifierPrepare <- function(counts, info, strand, segs, normVector = NULL, manual.segs=FALSE) {

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

  # check strand states
  assert_that(is.data.table(strand),
              "sample"%in% colnames(strand),
              "cell"  %in% colnames(strand),
              "chrom" %in% colnames(strand),
              "start" %in% colnames(strand),
              "end"   %in% colnames(strand),
              "class" %in% colnames(strand)) %>% invisible
  setkey(strand, sample, cell, chrom, start, end)

  # -> test whether strand states are disjoint intervals
  test_disjoint_intervals <- function(d) {
      gr <- makeGRangesFromDataFrame(d)
      end(gr) <- end(gr) -1
      return(isDisjoint(gr))
  }
  strand[, assert_that(test_disjoint_intervals(.SD)), by = .(sample, cell)] %>% invisible

  if (!manual.segs) {
    # segs
    assert_that(is.data.table(segs),
                "chrom" %in% colnames(segs),
                "bps"   %in% colnames(segs)) %>% invisible
    setkey(segs, chrom, bps)
  }


  # Kick out non-PASS cells
  if (nrow(info[pass1 != 1])> 0) {
    message("[MosaiClassifier] Kicking out ",
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

	if (manual.segs) {
		bin.size <- median(counts[,end-start])
    print(segs)
    print(info)
		probs <- merge(segs, info[,.(sample, cell, nb_p, mean)], by=c("sample", "cell"))
    print(probs)
    print(strand)
		message("[MosaiClassifier] Annotating strand-state")
    
		probs = addStrandStates(probs, strand)
		probs[, expected:=mean*(end-start)/bin.size, by=.(sample, cell)]
	}
	else {
		##############################################################################
		# Expand table to (all cells) x (all segments)
		#
		# Go from 0-based coordinates to 1-based bin coordinates
		segs[, to := bps + 1L]

		# add a "from" column which contians the "to" breakpoint from the prev. segment each
		segs[, from := (data.table::shift(to,fill = 0L) + 1L), by = chrom]

		# remove columns "bps" and "k"
		segs[, `:=`(bps = NULL, k = NULL)]




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
		message("[MosaiClassifier] Annotating observed W/C counts")
		probs <- addCountsPerSegment(probs, counts, manual.segs)
		probs[, `:=`(from = NULL,
               to   = NULL,)]
        }

  # Add normalization factors to the expected counts ("scalar")
  if (!is.null(normVector)) {
    message("[MosaiClassifier] Normalize coverage expectation")
    addNormalizationScalar(probs, counts, normVector)
  } else {
    probs[, scalar := 1.0]
  }


  probs[, mean := NULL]
  
  return(probs)
  
}




mosaiClassifierCalcProbs <- function(probs, maximumCN=4, haplotypeMode=F, alpha=0.05, definedHapStatus=FALSE, hapStatus=NULL, manual.segs=FALSE) {

  assert_that(is.data.table(probs),
              "sample" %in% colnames(probs),
              "cell"   %in% colnames(probs),
              "chrom"  %in% colnames(probs),
              "start"  %in% colnames(probs),
              "end"    %in% colnames(probs),
              "nb_p"   %in% colnames(probs),
              "expected"     %in% colnames(probs),
              "scalar" %in% colnames(probs),
              "W"      %in% colnames(probs),
              "C"      %in% colnames(probs),
              !("haplotype"  %in% colnames(probs)),
              !("haplo_name" %in% colnames(probs)),
              !("nb_hap_ll"  %in% colnames(probs)))

  if (!manual.segs){
    assert_that("num_bins" %in% colnames(probs))
  }


  if (!definedHapStatus){
    # defining the vector of all possible haplotypes
    hapStatus <- NULL
    for (j in 0:maximumCN){
      hapStatus <- c(hapStatus, allStatus(3,j))
    }
    for (j in 1:length(hapStatus)){
      hapStatus[j] <- paste(decodeStatus(hapStatus[j]), collapse = '')
    }
  }

  # creating a datatable containing all possible combinations of strand states and haplotypes,
  # and setting their segTypes
  hapStrandStates <- data.table()
  for (st in c("CC","WW","WC","CW")){
    hapStrandStates <- rbind(hapStrandStates, 
                             data.table(class=st, haplotype=hapStatus, 
                                        segtype=t(sapply(hapStatus, function(x) getSegType(st, x)))))
  }
  # naming third and forth columns
  colnames(hapStrandStates)[3:4]=c("Wcn", "Ccn")

  # Add SCE states where all possibilities are equal (by setting Wcn and Ccn to 0):
  hapStrandStates = rbind(hapStrandStates,
                          data.table(class = "?",
                                     haplotype = hapStatus,
                                     Wcn   = 1,
                                     Ccn   = 1))

  # Make sure that all states listed in probs are one of "WW","WC","CW","CC", or "?"
  assert_that(all(unique(probs$class) %in% c("WW","WC","CW","CC","?"))) %>% invisible

  # adding haplotype and genotype name columns
  hapStrandStates[,haplo_name:=.(sapply(haplotype, get_hap_name))]
  hapStrandStates[,geno_name:=.(sapply(haplo_name, haplo_to_geno_name))]

  # sort based on state
  setkey(hapStrandStates, class)

  ##### mering probs and haplotype strand states
  message("[MosaiClassifier] Expand table by ",
          length(hapStatus),
          " possible haplotype states")
  probs <- merge(probs, 
                 hapStrandStates,
                 by = "class",
                 allow.cartesian = T)

  ###########
  # reshuffling the columns
  if (manual.segs){
  probs <- probs[,.(sample, cell, chrom, start, end, class, nb_p, expected,
                    W, C, scalar, haplotype, Wcn, Ccn, haplo_name, geno_name)]
  } else {
  probs <- probs[,.(sample, cell, chrom, start, end, class, nb_p, expected, num_bins,
                    W, C, scalar, haplotype, Wcn, Ccn, haplo_name, geno_name)]
  }

  # computing dispersion parameters seperately for each segment and W and C counts ("Wcn" and "Ccn")
  message("[MosaiClassifier] Calculate dispersion parameters")

  probs <- add_dispPar(probs, alpha)

  # compute NB haplotype likelihoods
  message("[MosaiClassifier] Calculate NB likelihoods")
  probs[, nb_hap_ll := dnbinom(W, size = disp_w, prob = nb_p)
        * dnbinom(C, size = disp_c, prob = nb_p)]

  # computing sister haplotype (haplotype with the same genotype) for each haplotype
  sister.haps <- sapply(hapStatus, sisterHaplotype)
  sister.hap.pos <- match(sister.haps, hapStatus)
  # compute the set of symmetric haplotypes (haplotypes that are equal to their sister haplotype)
  symmetric.haps <- hapStatus[which(sister.haps==hapStatus)]

  # averaging the nb probs of sister haplotypes, when haplotype specific strand states are not known
  # adding genotype likelihoods, if haplotype mode is false
  if (!haplotypeMode) {
    message("[MosaiClassifier] Haplotype-unaware mode. H1 and H2 events are averaged.")
    probs[,nb_hap_ll:=.((nb_hap_ll+nb_hap_ll[sister.hap.pos])/2), by=.(sample, cell, chrom, start, end)]
  }

  # computing genotype likelihoods
  probs[,nb_gt_ll:=.(nb_hap_ll+nb_hap_ll[sister.hap.pos]), by=.(sample, cell, chrom, start, end)]
  # deviding the gt likelihoods of symmetric haplotypes by 2
  probs[haplotype %in% symmetric.haps, nb_gt_ll:=.(nb_gt_ll/2)]

}


# priors are 100 by default.
# Variant classes (matched by their `haplo_name`) which we want to weight differently
# can be included into the `priors` argument in form of a data table.
# A uniform prior can be done by specifying `priors = NULL`
mosaiClassifierPostProcessing <- function(probs,
                                          haplotypeMode=F,
                                          regularizationFactor=1e-10,
                                          priors = data.table(
                                                haplo_name = c("ref_hom","idup_h1","idup_h2","complex"),
                                                prior      = c(      200,       90,       90,        1)))
{
  assert_that(is.data.table(probs),
              "sample" %in% colnames(probs),
              "cell"   %in% colnames(probs),
              "chrom"  %in% colnames(probs),
              "start"  %in% colnames(probs),
              "end"    %in% colnames(probs),
              "haplo_name" %in% colnames(probs),
              "haplotype"  %in% colnames(probs),
              "nb_hap_ll"  %in% colnames(probs)) %>% invisible
  assert_that(!("nb_hap_pp" %in% colnames(probs))) %>% invisible

  message("[MosaiClassifier] regularizationFactor = ", regularizationFactor)

  # testing if there are some segments with zero probability for all haplotypes
  segs_max_hap_nb_probs <- probs[,
                                 .(sample, cell, chrom, start, end, max_nb_hap_ll=rep(max(nb_hap_ll), .N)),
                                 by=.(sample, cell, chrom, start, end)]
  message("[MosaiClassifier] The number of segments with 0 prob for all haplotypes is = ", 
          segs_max_hap_nb_probs[max_nb_hap_ll==0, .N])

  # set a uniform prob on sce segs and the segs_max_hap_nb_probs=0
  probs[segs_max_hap_nb_probs$max_nb_hap_ll==0,
        `:=`(nb_hap_ll = 1.0, nb_geno_ll = 1.0)]


  # PRIORS
  # when priors are given
  if (!is.null(priors)) {

    # first check that priors are in the right format
    assert_that(is.data.table(priors),
                "haplo_name" %in% colnames(priors),
                "prior"      %in% colnames(priors),
                all(is.numeric(priors$prior) & priors$prior >= 0),
                all(priors$haplo_name %in% probs$haplo_name)) %>% invisible

    # then add these priors to the probs table
    probs <- merge(probs, priors, by = "haplo_name", all.x = T)

    # set all other priors to 100
    probs[is.na(prior), prior := 100]

  # when no priors are given
  } else {
    message("[MosaiClassifier] Using uniform priors")
    probs[, prior:= 100]
  }


  # compute the posteriori probs (add new columns)
  probs[,nb_hap_pp:=.(nb_hap_ll*prior)][,nb_gt_pp:=.(nb_gt_ll*prior)]

  # normalizing nb_hap_pp and nb_gt_pp to 1 per sample, cell, and segment
  probs[, nb_hap_pp := nb_hap_pp/sum(nb_hap_pp), by=.(sample, cell, chrom, start, end)]
  probs[, nb_gt_pp := nb_gt_pp/sum(nb_gt_pp), by=.(sample, cell, chrom, start, end)]

  # get the number of haplotypes
  numHaps <- length(unique(probs$haplotype))
  # regularizing nb_hap_ll to set the min possible likelihood to a constant small number
  probs[,nb_hap_pp:=.((regularizationFactor/numHaps)+nb_hap_pp*(1-regularizationFactor))]
  probs[,nb_gt_pp:=.((regularizationFactor/numHaps)+nb_gt_pp*(1-regularizationFactor))]

  return(probs)
}
