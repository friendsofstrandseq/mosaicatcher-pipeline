library(dplyr)
library(data.table)
library(assertthat)
source("utils/mosaiClassifier/getStrandStates.R")
source("utils/mosaiClassifier/getCountsPerSegment.R")
source("utils/mosaiClassifier/generateHaploStates.R")
source("utils/mosaiClassifier/getDispParAndSegType.R")



getSVProbabilities <- function(counts, info, strand, segs) {

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
  probs <- add_seg_counts(probs, counts)
  probs[, scalar := 1]



  ##############################################################################
  # Extend table by haplotype states
  #
  @ work
  

  
