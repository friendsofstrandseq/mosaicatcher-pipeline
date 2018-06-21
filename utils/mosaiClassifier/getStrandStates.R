suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))

addPositions <- function(probs, counts) {
  assert_that(is.data.table(counts),
              "chrom" %in% colnames(counts),
              "start" %in% colnames(counts),
              "end"   %in% colnames(counts),
              "sample"%in% colnames(counts),
              "cell"  %in% colnames(counts)) %>% invisible
  assert_that(is.data.table(probs),
              "chrom" %in% colnames(probs),
              "from"  %in% colnames(probs),
              "to"    %in% colnames(probs))

  # get the list of bins
  bins <- unique(counts[, .(chrom_ = chrom, start, end)])
  setkey(bins, chrom_, start, end)

  # Make sure that all cells have exactly the same list of bins
  counts[,
         assert_that(all(unique(.SD) == bins)),
         by = .(sample,cell), .SDcols = c("chrom", "start", "end")] %>% invisible

  # Add start and end position
  probs[, `:=`(start = bins[chrom_ == chrom]$start[from],
               end   = bins[chrom_ == chrom]$end[to]),
        by = chrom]
  # Modified by reference (dummy return)
  probs
}


addStrandStates <- function(probs, strand) {

  assert_that(is.data.table(probs),
              "chrom" %in% colnames(probs),
              "start" %in% colnames(probs),
              "end"   %in% colnames(probs),
              "sample"%in% colnames(probs),
              "cell"  %in% colnames(probs)) %>% invisible
  assert_that(is.data.table(strand),
              "sample"%in% colnames(strand),
              "cell"  %in% colnames(strand),
              "chrom" %in% colnames(strand),
              "start" %in% colnames(strand),
              "end"   %in% colnames(strand),
              "class" %in% colnames(strand)) %>% invisible

  # Assign strand states to all chromosomes (ignoring the position, for now)
  xxx = merge(probs[, .(sample, cell, chrom, start, end)],
              strand[, .(sample, cell, chrom, class.start = start, class.end = end, class)],
              by = c("sample","cell","chrom"),
              all.x = T,
              allow.cartesian = T)

  # Keep only the segments/cells for which a strand state is annotated AND which
  # are fully contained within such a state region
  xxx <- xxx[!is.na(class) & start >= class.start & end <= class.end]

  # write annotation into the "probs" table for the cases which could successfully be annotated
  probs <- merge(probs,
                 xxx[, .(sample, cell, chrom, start, end, class)],
                 by = c("sample", "cell","chrom", "start", "end"),
                 all.x = T)

  # Set the class to "?" for the cases that have no correct annotation
  probs[is.na(class), class := "?"]
  message("[MosaiClassifier] Could assign strand states to ",
          nrow(xxx),
          " out of ",
          nrow(probs),
          " segments")
  # Return by reference
  probs
}
